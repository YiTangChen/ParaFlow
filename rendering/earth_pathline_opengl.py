#!/usr/bin/env python3
"""
Textured Earth in raw OpenGL — MPAS-Ocean scale (learning scaffold).

Draws a rotating, lit, texture-mapped globe using PyOpenGL for the GL calls
and GLFW for the window/context.  The sphere is built at the MPAS-Ocean
reference radius (6,371,229 m) so this can later host real MPAS pathline data
in the same coordinate system as rendering/render_3d_animation.py.

All seeds animate together, each leaving the growing pathline it has traced so
far, riding just above the globe surface.

Run:
    python rendering/opengl/earth_pathline_opengl.py
    python rendering/opengl/earth_pathline_opengl.py path/to/earth.jpg
    python rendering/opengl/earth_pathline_opengl.py path/to/earth.jpg path/to/pathlines.bin

Controls:
    drag left mouse = rotate      scroll = zoom        Space = toggle self-spin
    P = pause/resume              drag the bottom bar = scrub to any time
    left/right = step 1 hour      +Shift = 1 day       +Ctrl = 30 days
    G = go to a date (type e.g. 16-06-01 12:00, Enter to jump, Esc to cancel)
    Home/End = first/last step    R = restart          Esc = quit
    F = fade tails on/off         [ / ] = shorter/longer fade tails
    -/= = slower/faster           V = record MP4 on/off (needs ffmpeg on PATH)
An equirectangular "Blue Marble" JPEG gives a photoreal Earth; without one a
procedural land/ocean texture is generated so the program still runs.
"""

import sys
import os
import math
import struct
import ctypes
import shutil
import subprocess
from datetime import datetime, timedelta

import glfw
import numpy as np
from OpenGL.GL import *
from OpenGL.GLU import *

EARTH_RADIUS = 6_371_229.0   # metres — MPAS-Ocean reference sphere

# ---- Simulated clock ----
# pathlines.bin stores pure geometry — no timestamps — so the wall clock has to
# be reconstructed from the tracer settings: the integrator takes SOLVER_DT
# steps and write_trajectory keeps only every RECORD_INTERVAL-th one, so
# consecutive STORED points are SOLVER_DT * RECORD_INTERVAL apart.
SOLVER_DT       = 60.0                    # seconds per RK step (conf dt)
RECORD_INTERVAL = 60                      # steps between stored points (conf record_interval)
SAMPLE_SECONDS  = SOLVER_DT * RECORD_INTERVAL          # = 3600 s: one stored point per hour
SIM_START       = datetime(15, 1, 1)      # simulation year 0015, not calendar 2015
# The data window runs to 0017-12-31; traces end a little before that (the
# longest is 25,572 points = ~1,065 days), so dates stop short of the window's end.

# ---- Pathline animation tuning ----
DEFAULT_PATHLINES    = r"D:\Projects\pathline\highres_64_cpu_10k_2015_2017\pathlines.bin"
# Full resolution: every seed, every timestep.  For the 10k dataset that is
# 6,969 seeds / 81.4M points ~= 1.2 GB of GPU buffers, which is why the trails
# live in VBOs and are drawn with glMultiDrawArrays (see build_path_buffers).
# Raise either stride to trade fidelity back for memory.
SEED_STRIDE          = 1         # keep every Nth seed
TIME_STRIDE          = 1         # keep every Nth timestep along each path
LINE_LIFT            = 1.001     # draw trails at R*LIFT so they ride just above the globe
ANIM_STEPS_PER_FRAME = 40        # subsampled timesteps advanced per rendered frame
TRAIL_WIDTH          = 1.2
FADE_TIMESTEPS       = 1500      # real MPAS timesteps over which a trail fades out (time-based)
FADE_MIN_TIMESTEPS   = TIME_STRIDE * 2    # shortest fade window '[' can reach (2 samples)
FADE_MAX_TIMESTEPS   = 100_000            # longest fade window ']' can reach
# Fading tails are drawn as FADE_BANDS constant-alpha age bands rather than a
# per-vertex alpha array: the colours live in a static VBO, so the age ramp is
# applied with glBlendColor per band instead of re-uploading colours each frame.
FADE_BANDS           = 16

# ---- Timeline scrub bar (window pixels) ----
BAR_MARGIN = 90      # gap from the left/right window edges
BAR_BOTTOM = 26      # gap from the bottom window edge
BAR_HEIGHT = 10

# ---- Screen recording ('V') ----
REC_FPS = 60         # frames per second written into the MP4
REC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "recordings")
# Depth colour scale — trails are coloured by each point's depth below the
# surface (EARTH_RADIUS - |point|), interpolating SURFACE_COLOR -> DEEP_COLOR.
SURFACE_COLOR        = (1.0, 0.95, 0.15)   # yellow = sea surface (depth 0)
DEEP_COLOR           = (0.90, 0.05, 0.05)  # red    = at/below the colour-scale max
DEPTH_CLIM_PCT       = 95        # depth percentile of the data that saturates to red


# ---------------------------------------------------------------------------
# Geometry: a UV sphere in real metres
# ---------------------------------------------------------------------------
def make_sphere(radius, stacks=64, slices=64):
    """Return (verts Nx3, normals Nx3, uvs Nx2, indices M) for a UV sphere.

    Built vectorised with numpy.  Normals are the unit directions (position /
    radius) — free for a sphere and required by the lighting model.
    """
    i = np.arange(stacks + 1)
    j = np.arange(slices + 1)
    phi = (math.pi * i / stacks)[:, None]         # colatitude 0..pi (0 = north pole)
    theta = (2 * math.pi * j / slices)[None, :]   # longitude 0..2pi, east from +X

    # Pole along +Z and Greenwich (0 deg lon) along +X — the SAME Cartesian
    # convention MPAS uses (lat = asin(z/r), lon = atan2(y, x)).  Sharing this
    # frame is what makes the texture and the pathlines correspond exactly.
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi) * np.ones_like(theta)
    normals = np.stack([x, y, z], -1).reshape(-1, 3).astype(np.float32)
    verts = (normals * radius).astype(np.float32)

    # Blue-Marble maps are Greenwich-CENTRED (lon 0 at image middle, east to the
    # right), so u = lon/360 + 0.5.  Kept continuous (0.5..1.5) so GL_REPEAT
    # wraps the far side seamlessly.  v = 0 at the north pole (row 0 = North).
    u = (theta / (2 * math.pi) + 0.5) * np.ones_like(phi)
    v = (i / stacks)[:, None] * np.ones_like(theta)
    uvs = np.stack([u, v], -1).reshape(-1, 2).astype(np.float32)

    # Two triangles per grid quad, as indices into the vertex array.
    a = (i[:-1, None] * (slices + 1) + j[None, :-1])   # top-left of each quad
    b = a + (slices + 1)                               # one row down
    quad = np.stack([a, b, a + 1, a + 1, b, b + 1], -1)
    indices = quad.reshape(-1).astype(np.uint32)
    return verts, normals, uvs, indices


# ---------------------------------------------------------------------------
# Texture: load an equirectangular JPEG, or synthesize one
# ---------------------------------------------------------------------------
def procedural_earth(w=1024, h=512):
    """A crude land/ocean equirectangular image so the demo runs asset-free."""
    lat = np.linspace(90, -90, h)[:, None]        # rows: north -> south (matches v=0=N pole)
    lon = np.linspace(0, 360, w)[None, :]         # cols: 0 -> 360 east
    # Wavy pseudo-continents via summed sinusoids; > threshold = land.
    field = (np.sin(np.radians(lon * 3)) * np.cos(np.radians(lat * 2))
             + 0.6 * np.sin(np.radians(lon * 7 + lat * 4)))
    land = field > 0.35
    img = np.empty((h, w, 3), np.uint8)
    img[...] = (30, 90, 160)                       # ocean blue
    img[land] = (90, 130, 70)                      # land green
    img[np.abs(lat).repeat(w, 1) > 75] = (235, 235, 240)  # polar ice caps
    return img


def load_texture(path):
    """Upload an image to a GL texture. Returns the texture id."""
    if path:
        from PIL import Image
        pil = Image.open(path).convert("RGB")
        # Standard equirectangular maps have row 0 = North, which is exactly
        # what our v=0=north-pole convention wants — so DON'T flip vertically.
        img = np.asarray(pil, np.uint8)
    else:
        img = procedural_earth()          # also built with row 0 = North

    h, w = img.shape[:2]
    tex = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, tex)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0,
                 GL_RGB, GL_UNSIGNED_BYTE, np.ascontiguousarray(img))
    return tex


# ---------------------------------------------------------------------------
# Pathlines: load MPAS trajectories, subsample, project onto the globe
# ---------------------------------------------------------------------------
def _cache_path(path, seed_stride, time_stride):
    """Cache file for one (dataset, strides) pair, invalidated by size+mtime."""
    st = os.stat(path)
    return (f"{os.path.splitext(path)[0]}.s{seed_stride}_t{time_stride}"
            f"_{st.st_size}_{int(st.st_mtime)}.npz")


def read_pathlines(path, seed_stride, time_stride, radius):
    """Read pathlines.bin and return packed vertex data for animation.

    File format (little-endian): int32 n; then per pathline int32 pid,
    int32 nPts, nPts*3 float64 Cartesian metres.  The points already share the
    globe's coordinate system, so we only subsample and project each onto a
    sphere of `radius` (= EARTH_RADIUS * LINE_LIFT) so trails sit just above
    the surface instead of a few metres inside it.

    At full resolution the file is ~2 GB of float64, so it is parsed in two
    passes — headers first to size the output exactly, then one seed at a time
    straight into preallocated float32 buffers — rather than read whole and
    concatenated, which would need several GB of peak RAM.  The result is
    cached next to the dataset so only the first run pays for the parse.

    Returns (verts float32 (Ntot,3), offsets int32 (n,), lengths int32 (n,),
    max_len, rgb uint8 (Ntot,3), clim).  Seed i owns the slice
    verts[offsets[i] : offsets[i] + lengths[i]]; rgb aligns with verts and holds
    the depth colour (yellow at the surface -> red at/below `clim` metres).
    """
    cache = _cache_path(path, seed_stride, time_stride)
    if os.path.exists(cache):
        print("  reusing cache:", os.path.basename(cache))
        z = np.load(cache)
        lengths = z["lengths"]
        return (z["verts"], z["offsets"], lengths, int(lengths.max()),
                z["rgb"], float(z["clim"]))

    with open(path, "rb") as f:
        # Pass 1: walk the headers only (seeking over the point data) so we know
        # exactly how many points survive the strides before allocating.
        n_pl = struct.unpack("<i", f.read(4))[0]
        heads = []                                  # (file offset, nPts) per kept seed
        for k in range(n_pl):
            _pid, nPts = struct.unpack("<ii", f.read(8))
            if k % seed_stride == 0:
                heads.append((f.tell(), nPts))
            f.seek(nPts * 24, 1)

        lengths = np.array([-(-nPts // time_stride) for _, nPts in heads], np.int32)
        offsets = np.zeros(len(heads), np.int32)
        offsets[1:] = np.cumsum(lengths)[:-1]
        total = int(lengths.sum())
        print(f"  parsing {len(heads):,} seeds / {total:,} points "
              f"({total * 12 / 1e9:.2f} GB positions)...")

        # Pass 2: project each seed's points onto the lifted sphere in place.
        verts = np.empty((total, 3), np.float32)
        depths = np.empty(total, np.float32)
        for (doff, nPts), o, L in zip(heads, offsets.tolist(), lengths.tolist()):
            f.seek(doff)
            c = np.fromfile(f, "<f8", count=nPts * 3).reshape(nPts, 3)[::time_stride]
            r = np.sqrt(np.einsum("ij,ij->i", c, c))
            verts[o:o + L] = c / r[:, None] * radius
            depths[o:o + L] = np.maximum(EARTH_RADIUS - r, 0.0)   # metres below surface

    # Colour scale from a subsample — an exact percentile over 81M values would
    # sort a 300 MB copy for a number that only sets the legend's red point.
    clim = max(float(np.percentile(depths[::37], DEPTH_CLIM_PCT)), 1.0)
    surf = np.array(SURFACE_COLOR, np.float32) * 255.0
    deep = np.array(DEEP_COLOR, np.float32) * 255.0
    rgb = np.empty((total, 3), np.uint8)                # uint8: 1/4 the VRAM of float32
    for i in range(0, total, 4_000_000):                # chunked to bound temporaries
        d = np.clip(depths[i:i + 4_000_000] / clim, 0.0, 1.0)[:, None]
        rgb[i:i + 4_000_000] = ((1.0 - d) * surf + d * deep).astype(np.uint8)
    del depths

    try:
        np.savez(cache, verts=verts, offsets=offsets, lengths=lengths,
                 rgb=rgb, clim=clim)
        print("  cached to:", os.path.basename(cache))
    except OSError as e:
        print("  (cache write skipped:", e, ")")
    return verts, offsets, lengths, int(lengths.max()), rgb, clim


def sim_date(sample_idx):
    """Simulated date/time at a stored sample index, as 'YYYY-MM-DD HH:MM'.

    Formatted by hand rather than with strftime: '%Y' zero-pads inconsistently
    (and on Windows outright rejects) years below 1000, and these are year 15.
    """
    d = SIM_START + timedelta(seconds=sample_idx * TIME_STRIDE * SAMPLE_SECONDS)
    return (f"{d.year:04d}-{d.month:02d}-{d.day:02d} "
            f"{d.hour:02d}:{d.minute:02d}")


def parse_sim_date(text):
    """'YYYY-MM-DD [HH[:MM]]' -> stored-sample index, or None if unparseable.

    The inverse of sim_date().  Years are simulation years (15..17), so '15-6-1'
    and '0015-06-01 12:30' both work.  The returned index may fall outside the
    data — the caller clamps.
    """
    parts = text.strip().replace("/", "-").split()
    if not parts:
        return None
    ymd = parts[0].split("-")
    if len(ymd) != 3:
        return None
    try:
        hh, mm = 0, 0
        if len(parts) > 1:
            hms = parts[1].split(":")
            hh = int(hms[0])
            mm = int(hms[1]) if len(hms) > 1 else 0
        when = datetime(int(ymd[0]), int(ymd[1]), int(ymd[2]), hh, mm)
    except ValueError:
        return None                      # bad numbers, month 13, 31 Feb, ...
    secs = (when - SIM_START).total_seconds()
    return secs / (TIME_STRIDE * SAMPLE_SECONDS)


def seek_to(sample_idx):
    """Jump the animation to a stored-sample index, clamped to the data."""
    if View.max_len < 1:
        return
    idx = max(0, min(View.max_len - 1, int(sample_idx)))
    View.t_accum = float(idx)            # keep the fractional clock in step
    View.t_anim = idx


def timeline_rect(win):
    """Scrub-bar rectangle in WINDOW coords (x0, y0, x1, y1), y from the top.

    Window (not framebuffer) coords, because that is what glfw.get_cursor_pos
    reports — on a HiDPI display the two differ by the content scale.
    """
    w, h = glfw.get_window_size(win)
    return BAR_MARGIN, h - BAR_BOTTOM - BAR_HEIGHT, w - BAR_MARGIN, h - BAR_BOTTOM


def draw_timeline(win, t):
    """Draw the scrub bar: track, elapsed fill, and a knob at the current time."""
    x0, y0, x1, y1 = timeline_rect(win)
    ww, wh = glfw.get_window_size(win)
    fbw, fbh = glfw.get_framebuffer_size(win)
    s = fbw / max(ww, 1)                              # window px -> framebuffer px
    # Flip y (window y counts down from the top, GL counts up from the bottom).
    bx0, bx1 = x0 * s, x1 * s
    by0, by1 = fbh - y1 * s, fbh - y0 * s
    frac = t / max(View.max_len - 1, 1)
    knob = bx0 + (bx1 - bx0) * frac

    glMatrixMode(GL_PROJECTION)
    glPushMatrix()
    glLoadIdentity()
    glOrtho(0, fbw, 0, fbh, -1, 1)                    # pixel-space overlay
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    glDisable(GL_DEPTH_TEST)
    glDisable(GL_TEXTURE_2D)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def quad(qx0, qy0, qx1, qy1, rgba):
        glColor4f(*rgba)
        glBegin(GL_QUADS)
        glVertex2f(qx0, qy0); glVertex2f(qx1, qy0)
        glVertex2f(qx1, qy1); glVertex2f(qx0, qy1)
        glEnd()

    quad(bx0, by0, bx1, by1, (0.0, 0.0, 0.0, 0.45))              # track
    quad(bx0, by0, knob, by1, (1.0, 0.92, 0.35, 0.85))           # elapsed
    kw = 3.0 * s
    quad(knob - kw, by0 - 4 * s, knob + kw, by1 + 4 * s,
         (1.0, 1.0, 1.0, 0.95))                                  # knob

    glDisable(GL_BLEND)
    glEnable(GL_DEPTH_TEST)
    glColor3f(1, 1, 1)                                # leave colour state clean
    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glPopMatrix()
    glMatrixMode(GL_MODELVIEW)


def build_path_buffers(verts, rgb):
    """Upload trail positions/colours to VBOs once; return (pos_vbo, col_vbo).

    At full resolution these are ~1.0 GB of positions and ~0.24 GB of colours.
    Client-side arrays would re-stream all of that across the bus every frame,
    so the data lives in GPU memory and each frame only issues draw ranges.
    """
    pos_vbo, col_vbo = glGenBuffers(2)
    glBindBuffer(GL_ARRAY_BUFFER, pos_vbo)
    glBufferData(GL_ARRAY_BUFFER, verts.nbytes, verts, GL_STATIC_DRAW)
    glBindBuffer(GL_ARRAY_BUFFER, col_vbo)
    glBufferData(GL_ARRAY_BUFFER, rgb.nbytes, rgb, GL_STATIC_DRAW)
    glBindBuffer(GL_ARRAY_BUFFER, 0)
    err = glGetError()
    if err == GL_OUT_OF_MEMORY:
        raise MemoryError(
            f"GPU could not hold {(verts.nbytes + rgb.nbytes) / 1e9:.2f} GB of trails — "
            f"raise SEED_STRIDE or TIME_STRIDE")
    print(f"  uploaded {(verts.nbytes + rgb.nbytes) / 1e9:.2f} GB to GPU buffers")
    return pos_vbo, col_vbo


def draw_line_strips(firsts, counts):
    """One glMultiDrawArrays over every seed's visible range (skipping empties)."""
    keep = counts >= 2                       # a 1-point strip draws nothing
    if not keep.any():
        return
    f = np.ascontiguousarray(firsts[keep], np.int32)
    c = np.ascontiguousarray(counts[keep], np.int32)
    glMultiDrawArrays(GL_LINE_STRIP, f, c, len(f))


# ---------------------------------------------------------------------------
# Text overlay — rasterised with Pillow, blitted with glDrawPixels
# ---------------------------------------------------------------------------
def draw_text_overlay(win, text, cache, top=12, fill=(255, 235, 90, 255)):
    """Blit `text` at the top-left in window space, `top` px down from the top.

    `cache` (a dict) persists the font and last-rendered image so we only
    re-rasterise when the text actually changes.  glDrawPixels bypasses the
    modelview/projection matrices, so this is independent of the 3D camera.
    """
    if cache.get("text") != text:
        from PIL import Image, ImageDraw, ImageFont
        font = cache.get("font")
        if font is None:
            for name in ("arialbd.ttf", "arial.ttf", "DejaVuSans-Bold.ttf"):
                try:
                    font = ImageFont.truetype(name, 22)
                    break
                except OSError:
                    continue
            cache["font"] = font or ImageFont.load_default()
            font = cache["font"]
        l, t, r, b = ImageDraw.Draw(Image.new("RGBA", (1, 1))).textbbox((0, 0), text, font=font)
        w, h = (r - l) + 12, (b - t) + 10
        img = Image.new("RGBA", (w, h), (0, 0, 0, 110))        # translucent dark plate
        ImageDraw.Draw(img).text((6 - l, 5 - t), text, font=font, fill=fill)
        # glDrawPixels draws the bottom row first, so flip vertically.
        cache["arr"] = np.ascontiguousarray(np.asarray(img, np.uint8)[::-1])
        cache["w"], cache["h"], cache["text"] = w, h, text

    arr, w, h = cache["arr"], cache["w"], cache["h"]
    _, fbh = glfw.get_framebuffer_size(win)
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glWindowPos2i(12, fbh - h - top)                           # top-left (y from bottom)
    glDrawPixels(w, h, GL_RGBA, GL_UNSIGNED_BYTE, arr)
    glDisable(GL_BLEND)
    glEnable(GL_DEPTH_TEST)


def draw_depth_legend(win, clim, cache):
    """Blit a vertical yellow->red depth colourbar on the right edge.

    Built once (into `cache`) since the scale is fixed after load: a translucent
    plate with a title, the gradient bar (surface at top, deep at bottom) and
    tick labels at 0, clim/2 and clim.
    """
    if "arr" not in cache:
        from PIL import Image, ImageDraw, ImageFont
        font = None
        for name in ("arial.ttf", "DejaVuSans.ttf"):
            try:
                font = ImageFont.truetype(name, 16)
                break
            except OSError:
                continue
        font = font or ImageFont.load_default()

        bar_w, bar_h, pad_top = 22, 200, 26
        w, h = 12 + bar_w + 10 + 56, pad_top + bar_h + 14
        img = Image.new("RGBA", (w, h), (0, 0, 0, 120))         # translucent plate
        d = ImageDraw.Draw(img)
        d.text((10, 5), "Depth (m)", font=font, fill=(255, 255, 255, 255))
        x0, y0 = 12, pad_top
        surf, deep = np.array(SURFACE_COLOR), np.array(DEEP_COLOR)
        for i in range(bar_h):                                  # surface (top) -> deep (bottom)
            c = surf * (1 - i / (bar_h - 1)) + deep * (i / (bar_h - 1))
            d.line([(x0, y0 + i), (x0 + bar_w, y0 + i)],
                   fill=(int(c[0] * 255), int(c[1] * 255), int(c[2] * 255), 255))
        d.rectangle([x0, y0, x0 + bar_w, y0 + bar_h], outline=(230, 230, 230, 220))
        for frac, lab in ((0.0, "0"), (0.5, f"{clim * 0.5:.0f}"), (1.0, f"{clim:.0f}+")):
            yy = y0 + int(frac * (bar_h - 1))
            d.line([(x0 + bar_w, yy), (x0 + bar_w + 4, yy)], fill=(230, 230, 230, 220))
            d.text((x0 + bar_w + 8, yy - 8), lab, font=font, fill=(255, 255, 255, 255))
        cache["arr"] = np.ascontiguousarray(np.asarray(img, np.uint8)[::-1])
        cache["w"], cache["h"] = w, h

    arr, w, h = cache["arr"], cache["w"], cache["h"]
    fbw, fbh = glfw.get_framebuffer_size(win)
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glWindowPos2i(fbw - w - 14, (fbh - h) // 2)                 # right edge, vertically centred
    glDrawPixels(w, h, GL_RGBA, GL_UNSIGNED_BYTE, arr)
    glDisable(GL_BLEND)
    glEnable(GL_DEPTH_TEST)


# ---------------------------------------------------------------------------
# Screen recorder: pipe the rendered back buffer straight to an MP4 via ffmpeg
# (or dump PNG frames if ffmpeg isn't on PATH).  One rendered frame == one video
# frame, so the clip plays at REC_FPS regardless of the live render rate.
# ---------------------------------------------------------------------------
class Recorder:
    def __init__(self):
        self.active = False
        self.proc = None            # ffmpeg subprocess, or None in PNG mode
        self.png_dir = None
        self.path = None
        self.w = self.h = 0
        self.frames = 0

    def toggle(self, win):
        self.stop() if self.active else self.start(win)

    def start(self, win):
        w, h = glfw.get_framebuffer_size(win)
        w -= w & 1                  # yuv420p needs even dimensions
        h -= h & 1
        if w < 2 or h < 2:
            return
        os.makedirs(REC_DIR, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.w, self.h, self.frames = w, h, 0
        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg:
            self.path = os.path.join(REC_DIR, f"mpaso_{stamp}.mp4")
            self.proc = subprocess.Popen(
                [ffmpeg, "-y", "-f", "rawvideo", "-pixel_format", "rgb24",
                 "-video_size", f"{w}x{h}", "-framerate", str(REC_FPS), "-i", "-",
                 "-an", "-c:v", "libx264", "-preset", "fast", "-crf", "18",
                 "-pix_fmt", "yuv420p", self.path],
                stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.png_dir = None
            print(f"● recording {w}x{h} @ {REC_FPS}fps -> {self.path}")
        else:
            self.proc = None
            self.png_dir = os.path.join(REC_DIR, f"mpaso_{stamp}_frames")
            os.makedirs(self.png_dir, exist_ok=True)
            print(f"● ffmpeg not found; recording PNG frames -> {self.png_dir}")
        self.active = True

    def capture(self):
        """Read the just-drawn back buffer and feed one frame to the encoder."""
        if not self.active:
            return
        glReadBuffer(GL_BACK)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        data = glReadPixels(0, 0, self.w, self.h, GL_RGB, GL_UNSIGNED_BYTE)
        arr = (np.frombuffer(data, np.uint8) if isinstance(data, (bytes, bytearray))
               else np.asarray(data, np.uint8)).reshape(self.h, self.w, 3)
        arr = arr[::-1]             # GL origin is bottom-left; flip to top-down
        if self.proc is not None:
            try:
                self.proc.stdin.write(np.ascontiguousarray(arr).tobytes())
            except (BrokenPipeError, OSError):
                print("  recording stopped: ffmpeg pipe closed")
                self.active = False
                self.proc = None
                return
        else:
            from PIL import Image
            Image.fromarray(arr).save(
                os.path.join(self.png_dir, f"frame_{self.frames:05d}.png"))
        self.frames += 1

    def stop(self):
        if not self.active:
            return
        self.active = False
        if self.proc is not None:
            try:
                self.proc.stdin.close()
            except OSError:
                pass
            self.proc.wait()
            print(f"■ saved {self.frames} frames -> {self.path}")
            self.proc = None
        else:
            print(f"■ saved {self.frames} PNG frames -> {self.png_dir}")


REC = Recorder()


# ---------------------------------------------------------------------------
# Interaction state (module-level so GLFW callbacks can mutate it)
# ---------------------------------------------------------------------------
class View:
    yaw = 0.0                    # degrees, mouse-controlled spin
    pitch = 0.0                  # degrees, tilt
    dist = EARTH_RADIUS * 2.8    # camera distance in metres
    dragging = False
    last = (0.0, 0.0)
    auto_rotate = False          # idle self-spin, toggled with Space
    t_anim = 0                   # current animation timestep (subsampled)
    t_accum = 0.0                # fractional accumulator so speed can be < 1 step/frame
    speed = float(ANIM_STEPS_PER_FRAME)   # subsampled steps advanced per frame ('-' / '=')
    fade = False                 # 'F' toggles fading comet tails vs persistent trails
    fade_steps = float(FADE_TIMESTEPS)    # tail length in real timesteps ('[' / ']')
    paused = False               # 'P' / Space-bar of playback: freeze the clock
    max_len = 0                  # samples in the longest path; set once loaded
    scrubbing = False            # dragging the timeline bar
    typing = False               # 'G' opens a date-entry line
    type_buf = ""                # what has been typed so far
    type_msg = ""                # validation feedback for the entry line


def scrub_at(win, x):
    """Seek to the time under cursor x (window coords) on the scrub bar."""
    x0, _, x1, _ = timeline_rect(win)
    frac = (x - x0) / max(x1 - x0, 1)
    seek_to(round(max(0.0, min(1.0, frac)) * (View.max_len - 1)))


def on_mouse_button(win, button, action, mods):
    if button != glfw.MOUSE_BUTTON_LEFT:
        return
    if action == glfw.PRESS:
        x, y = glfw.get_cursor_pos(win)
        x0, y0, x1, y1 = timeline_rect(win)
        # Generous vertical grab margin — the bar itself is only 10 px tall.
        if View.max_len > 1 and x0 - 8 <= x <= x1 + 8 and y0 - 8 <= y <= y1 + 8:
            View.scrubbing = True        # grabbed the timeline, not the globe
            scrub_at(win, x)
            return
        View.dragging = True
        View.last = (x, y)
    else:
        View.dragging = False
        View.scrubbing = False


def on_cursor(win, x, y):
    if View.scrubbing:
        scrub_at(win, x)
    elif View.dragging:
        dx, dy = x - View.last[0], y - View.last[1]
        View.yaw += dx * 0.3
        View.pitch = max(-89, min(89, View.pitch + dy * 0.3))
        View.last = (x, y)


def on_char(win, codepoint):
    """Collect date-entry keystrokes (digits and separators only)."""
    if View.typing and len(View.type_buf) < 20:
        ch = chr(codepoint)
        if ch in "0123456789-: /":
            View.type_buf += ch
            View.type_msg = ""


def on_scroll(win, xoff, yoff):
    View.dist *= 0.9 ** yoff      # scroll up = zoom in
    # Never let the camera reach the surface (radius R) — clamp just above it.
    View.dist = max(EARTH_RADIUS * 1.05, min(EARTH_RADIUS * 8, View.dist))


def on_key(win, key, scancode, action, mods):
    if action not in (glfw.PRESS, glfw.REPEAT):    # REPEAT: holding an arrow scrubs
        return

    # --- date entry swallows every key while open ---
    if View.typing:
        if key == glfw.KEY_ESCAPE:
            View.typing, View.type_buf, View.type_msg = False, "", ""
        elif key == glfw.KEY_BACKSPACE:
            View.type_buf = View.type_buf[:-1]
        elif key in (glfw.KEY_ENTER, glfw.KEY_KP_ENTER):
            idx = parse_sim_date(View.type_buf)
            if idx is None:
                View.type_msg = "  <- can't read that"
            elif not (0 <= idx <= View.max_len - 1):
                View.type_msg = f"  <- outside {sim_date(0)} .. {sim_date(View.max_len - 1)}"
            else:
                seek_to(idx)
                View.paused = True                 # land on the date and hold there
                View.typing, View.type_buf, View.type_msg = False, "", ""
        return

    if key == glfw.KEY_ESCAPE:
        glfw.set_window_should_close(win, True)
    elif key == glfw.KEY_SPACE:
        View.auto_rotate = not View.auto_rotate    # start/stop the self-rotation
    elif key == glfw.KEY_P:
        View.paused = not View.paused              # stop/resume the clock
    elif key == glfw.KEY_G:
        View.typing, View.type_buf, View.type_msg = True, "", ""   # go to a date
    elif key == glfw.KEY_V:
        REC.toggle(win)                            # start/stop MP4 recording
    elif key in (glfw.KEY_LEFT, glfw.KEY_RIGHT):
        # Step the clock: 1 hour, Shift = 1 day, Ctrl = 30 days.
        step = 1
        if mods & glfw.MOD_SHIFT:
            step = int(86400 / (TIME_STRIDE * SAMPLE_SECONDS)) or 1
        if mods & glfw.MOD_CONTROL:
            step = int(30 * 86400 / (TIME_STRIDE * SAMPLE_SECONDS)) or 1
        seek_to(View.t_anim + (step if key == glfw.KEY_RIGHT else -step))
        View.paused = True                         # stepping implies you want it still
    elif key == glfw.KEY_HOME:
        seek_to(0)
    elif key == glfw.KEY_END:
        seek_to(View.max_len - 1)
    elif key == glfw.KEY_R:
        View.t_anim = 0                            # restart the pathline animation
        View.t_accum = 0.0
    elif key == glfw.KEY_F:
        View.fade = not View.fade                  # fading tails <-> persistent trails
    elif key in (glfw.KEY_MINUS, glfw.KEY_COMMA):
        View.speed = max(0.25, View.speed * 0.5)   # slow down
    elif key in (glfw.KEY_EQUAL, glfw.KEY_PERIOD):
        View.speed = min(512.0, View.speed * 2.0)  # speed up (25,572 steps at full res)
    elif key == glfw.KEY_LEFT_BRACKET:             # shorter fading tails
        View.fade_steps = max(FADE_MIN_TIMESTEPS, View.fade_steps * 0.5)
    elif key == glfw.KEY_RIGHT_BRACKET:            # longer fading tails
        View.fade_steps = min(FADE_MAX_TIMESTEPS, View.fade_steps * 2.0)


# ---------------------------------------------------------------------------
def main():
    tex_path = sys.argv[1] if len(sys.argv) > 1 else None

    if not glfw.init():
        raise RuntimeError("GLFW failed to initialise")
    win = glfw.create_window(960, 720, "MPAS-O Earth (OpenGL)", None, None)
    if not win:
        glfw.terminate()
        raise RuntimeError("Failed to create GLFW window")
    glfw.make_context_current(win)
    glfw.swap_interval(1)          # vsync
    glfw.set_mouse_button_callback(win, on_mouse_button)
    glfw.set_cursor_pos_callback(win, on_cursor)
    glfw.set_scroll_callback(win, on_scroll)
    glfw.set_key_callback(win, on_key)
    glfw.set_char_callback(win, on_char)

    glEnable(GL_DEPTH_TEST)
    glEnable(GL_TEXTURE_2D)
    # No lighting: the texture is drawn at full, uniform brightness everywhere
    # (constant illumination, no day/night terminator).

    verts, normals, uvs, indices = make_sphere(EARTH_RADIUS)
    texture = load_texture(tex_path)
    print("Texture:", tex_path if tex_path else
          "(procedural — pass a Blue Marble JPEG for photoreal)")

    # --- Load pathlines (optional) ---
    pl_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_PATHLINES
    offsets = lengths = None
    pos_vbo = col_vbo = None
    max_len = 0
    clim = 0.0
    if pl_path and os.path.exists(pl_path):
        print("Loading pathlines:", pl_path)
        path_verts, offsets, lengths, max_len, rgb, clim = read_pathlines(
            pl_path, SEED_STRIDE, TIME_STRIDE, EARTH_RADIUS * LINE_LIFT)
        print(f"  {len(lengths):,} seeds (every {SEED_STRIDE}), "
              f"{path_verts.shape[0]:,} points, {max_len:,} animation steps")
        print(f"  depth colour saturates at p{DEPTH_CLIM_PCT} = {clim:.0f} m")
        print(f"  clock: 1 stored point = {SOLVER_DT:g}s x {RECORD_INTERVAL} steps "
              f"= {SAMPLE_SECONDS / 3600:g} h")
        print(f"  spans {sim_date(0)} -> {sim_date(max_len - 1)} "
              f"({(max_len - 1) * TIME_STRIDE * SAMPLE_SECONDS / 86400:.1f} days)")
        pos_vbo, col_vbo = build_path_buffers(path_verts, rgb)
        del path_verts, rgb            # the GPU owns the only copy from here on
        View.max_len = max_len         # callbacks need it to clamp seeks
    else:
        print("No pathlines file — showing globe only. (looked for", pl_path, ")")

    glEnableClientState(GL_VERTEX_ARRAY)

    spin = 0.0
    text_cache = {}
    legend_cache = {}
    entry_cache = {}

    while not glfw.window_should_close(win):
        # Fade window in TIME: how many stored samples span the current tail
        # length (View.fade_steps, adjustable with '[' / ']').
        fade_pts = max(2, int(View.fade_steps) // TIME_STRIDE)

        fbw, fbh = glfw.get_framebuffer_size(win)
        glViewport(0, 0, fbw, fbh)
        glClearColor(0.02, 0.02, 0.05, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # --- projection: near/far scaled to metres to keep depth buffer sane ---
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        # Near/far track the camera distance so the near plane sits just in
        # front of the globe's nearest point (dist - R) and never clips into it.
        near = max((View.dist - EARTH_RADIUS) * 0.5, EARTH_RADIUS * 1e-3)
        far = View.dist + EARTH_RADIUS * 2.0
        gluPerspective(45.0, fbw / max(fbh, 1), near, far)

        # --- camera / model transform ---
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -View.dist)
        glRotatef(View.pitch, 1, 0, 0)
        glRotatef(View.yaw,   0, 1, 0)
        glRotatef(-90, 1, 0, 0)              # +Z pole (MPAS north) points up
        glRotatef(spin, 0, 0, 1)             # self-spin about the true pole (+Z)

        # --- draw the textured globe (client arrays; one glDrawElements) ---
        # With a VBO bound, gl*Pointer would read its argument as a buffer
        # offset instead of a pointer, so unbind before the globe's arrays.
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_TEXTURE_COORD_ARRAY)
        glEnable(GL_TEXTURE_2D)
        glVertexPointer(3, GL_FLOAT, 0, verts)
        glNormalPointer(GL_FLOAT, 0, normals)
        glTexCoordPointer(2, GL_FLOAT, 0, uvs)
        glBindTexture(GL_TEXTURE_2D, texture)
        glColor3f(1, 1, 1)
        glDrawElements(GL_TRIANGLES, len(indices), GL_UNSIGNED_INT, indices)

        # --- draw pathlines: growing trails (all from GPU-resident VBOs) ---
        if pos_vbo is not None:
            # Only VERTEX_ARRAY/COLOR_ARRAY should be active, else GL reads
            # normals/uvs (sized for the globe) out of bounds for the trails.
            glDisableClientState(GL_NORMAL_ARRAY)
            glDisableClientState(GL_TEXTURE_COORD_ARRAY)
            glDisable(GL_TEXTURE_2D)
            glEnableClientState(GL_COLOR_ARRAY)
            glBindBuffer(GL_ARRAY_BUFFER, pos_vbo)
            glVertexPointer(3, GL_FLOAT, 0, ctypes.c_void_p(0))
            glBindBuffer(GL_ARRAY_BUFFER, col_vbo)
            glColorPointer(3, GL_UNSIGNED_BYTE, 0, ctypes.c_void_p(0))

            t = View.t_anim
            glLineWidth(TRAIL_WIDTH)
            newest = np.minimum(t, lengths - 1)      # newest existing sample per seed

            if View.fade:
                # Time-based fade: a point's alpha follows its age (t - step) in the
                # GLOBAL timeline, so even finished seeds keep fading and vanish once
                # older than the fade window.  The colours are static in a VBO, so the
                # ramp is quantised into FADE_BANDS constant-alpha slabs and applied
                # with glBlendColor — no per-frame colour uploads.
                glEnable(GL_BLEND)
                glBlendFunc(GL_CONSTANT_ALPHA, GL_ONE_MINUS_CONSTANT_ALPHA)
                fp1 = fade_pts - 1
                for band in range(FADE_BANDS):
                    age_hi = fp1 * band / FADE_BANDS          # age of the band's newest point
                    age_lo = fp1 * (band + 1) / FADE_BANDS    # age of the band's oldest point
                    alpha = 1.0 - (age_hi + age_lo) * 0.5 / fp1
                    if alpha <= 0.02:
                        continue                              # already invisible
                    # Band = global steps [t-age_lo, t-age_hi], clipped per seed. The
                    # +1 makes neighbouring bands share a vertex so the line has no gaps.
                    start = np.maximum(t - int(age_lo), 0)
                    end = np.minimum(newest, t - int(age_hi) + 1)
                    glBlendColor(0.0, 0.0, 0.0, alpha)
                    draw_line_strips(offsets + start, end - start + 1)
                glDisable(GL_BLEND)
            else:
                # Persistent trails: the whole path 0..t, coloured by depth.
                draw_line_strips(offsets, newest + 1)

            glBindBuffer(GL_ARRAY_BUFFER, 0)
            glDisableClientState(GL_COLOR_ARRAY)

            # Simulated clock, top-left, with the raw step for cross-checking.
            real_t = View.t_anim * TIME_STRIDE
            real_max = (max_len - 1) * TIME_STRIDE
            tail_hours = fade_pts * TIME_STRIDE * SAMPLE_SECONDS / 3600.0
            mode = (f"fading ({tail_hours / 24:.1f} d)"
                    if View.fade else "full")
            state = "PAUSED" if View.paused else f"speed {View.speed:g}"
            if REC.active:
                state += f"    REC {REC.frames}"
            draw_text_overlay(win,
                              f"{sim_date(View.t_anim)}    "
                              f"step {real_t:,} / {real_max:,}    trails: {mode}"
                              f"    {state}",
                              text_cache)
            if View.typing:
                draw_text_overlay(win,
                                  f"go to date: {View.type_buf}_"
                                  f"{View.type_msg or '   (YYYY-MM-DD [HH:MM], Enter/Esc)'}",
                                  entry_cache, top=54, fill=(140, 255, 170, 255))
            draw_depth_legend(win, clim, legend_cache)
            draw_timeline(win, View.t_anim)

        REC.capture()             # grab the finished back buffer (incl. HUD) if recording
        glfw.swap_buffers(win)
        glfw.poll_events()
        if View.auto_rotate and not View.dragging:
            spin += 0.15          # eastward (real Earth spin): features drift right
        # Playback only advances when running — pausing, scrubbing and date entry
        # all hold the clock wherever the user put it.
        if (pos_vbo is not None and max_len > 1
                and not (View.paused or View.scrubbing or View.typing)):
            View.t_accum += View.speed          # fractional advance (supports slow speeds)
            if View.t_accum >= max_len:
                View.t_accum = 0.0              # loop the animation
            View.t_anim = int(View.t_accum)

    REC.stop()                    # finalise the MP4 if a recording is still running
    glfw.terminate()


if __name__ == "__main__":
    main()
