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
    drag left mouse = rotate   scroll = zoom   Space = toggle self-spin
    F = fade tails on/off      -/= = slower/faster   R = restart   Esc = quit
An equirectangular "Blue Marble" JPEG gives a photoreal Earth; without one a
procedural land/ocean texture is generated so the program still runs.
"""

import sys
import os
import math
import struct

import glfw
import numpy as np
from OpenGL.GL import *
from OpenGL.GLU import *

EARTH_RADIUS = 6_371_229.0   # metres — MPAS-Ocean reference sphere

# ---- Pathline animation tuning ----
DEFAULT_PATHLINES    = r"D:\Projects\pathline\nersc_highres_16_cpu_10k\pathlines.bin"
SEED_STRIDE          = 10        # keep every Nth seed (6969 -> ~700) for memory/speed
TIME_STRIDE          = 10        # keep every Nth timestep along each path
LINE_LIFT            = 1.001     # draw trails at R*LIFT so they ride just above the globe
ANIM_STEPS_PER_FRAME = 4         # subsampled timesteps advanced per rendered frame
TRAIL_WIDTH          = 1.2
FADE_TIMESTEPS       = 1500      # real MPAS timesteps over which a trail fades out (time-based)
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
def read_pathlines(path, seed_stride, time_stride, radius):
    """Read pathlines.bin and return packed vertex data for animation.

    File format (little-endian): int32 n; then per pathline int32 pid,
    int32 nPts, nPts*3 float64 Cartesian metres.  The points already share the
    globe's coordinate system, so we only subsample and project each onto a
    sphere of `radius` (= EARTH_RADIUS * LINE_LIFT) so trails sit just above
    the surface instead of a few metres inside it.

    Returns (verts float32 (Ntot,3), offsets int32 (n,), lengths int32 (n,),
    max_len, depths float32 (Ntot,)).  Seed i owns the slice
    verts[offsets[i] : offsets[i] + lengths[i]]; depths align with verts and
    hold each point's metres below the surface (before the surface projection).
    """
    with open(path, "rb") as f:
        data = f.read()
    off = 0
    n_pl = struct.unpack_from("<i", data, off)[0]
    off += 4
    segs = []
    for k in range(n_pl):
        _pid, nPts = struct.unpack_from("<ii", data, off)
        off += 8
        if k % seed_stride == 0:
            c = np.frombuffer(data, "<f8", count=nPts * 3, offset=off).reshape(nPts, 3)
            segs.append(c[::time_stride])
        off += nPts * 24

    lengths = np.array([len(s) for s in segs], dtype=np.int32)
    offsets = np.zeros(len(segs), dtype=np.int32)
    offsets[1:] = np.cumsum(lengths)[:-1]

    allpts = np.concatenate(segs, axis=0)
    r = np.sqrt(np.einsum("ij,ij->i", allpts, allpts))
    depths = np.maximum(EARTH_RADIUS - r, 0.0).astype(np.float32)   # metres below surface
    verts = (allpts / r[:, None] * radius).astype(np.float32)
    return verts, offsets, lengths, int(lengths.max()), depths


# ---------------------------------------------------------------------------
# Text overlay — rasterised with Pillow, blitted with glDrawPixels
# ---------------------------------------------------------------------------
def draw_text_overlay(win, text, cache):
    """Blit `text` at the top-left in window space.

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
        ImageDraw.Draw(img).text((6 - l, 5 - t), text, font=font, fill=(255, 235, 90, 255))
        # glDrawPixels draws the bottom row first, so flip vertically.
        cache["arr"] = np.ascontiguousarray(np.asarray(img, np.uint8)[::-1])
        cache["w"], cache["h"], cache["text"] = w, h, text

    arr, w, h = cache["arr"], cache["w"], cache["h"]
    _, fbh = glfw.get_framebuffer_size(win)
    glDisable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    glWindowPos2i(12, fbh - h - 12)                            # top-left (y from bottom)
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


def on_mouse_button(win, button, action, mods):
    if button == glfw.MOUSE_BUTTON_LEFT:
        View.dragging = (action == glfw.PRESS)
        View.last = glfw.get_cursor_pos(win)


def on_cursor(win, x, y):
    if View.dragging:
        dx, dy = x - View.last[0], y - View.last[1]
        View.yaw += dx * 0.3
        View.pitch = max(-89, min(89, View.pitch + dy * 0.3))
        View.last = (x, y)


def on_scroll(win, xoff, yoff):
    View.dist *= 0.9 ** yoff      # scroll up = zoom in
    # Never let the camera reach the surface (radius R) — clamp just above it.
    View.dist = max(EARTH_RADIUS * 1.05, min(EARTH_RADIUS * 8, View.dist))


def on_key(win, key, scancode, action, mods):
    if action != glfw.PRESS:
        return
    if key == glfw.KEY_ESCAPE:
        glfw.set_window_should_close(win, True)
    elif key == glfw.KEY_SPACE:
        View.auto_rotate = not View.auto_rotate    # start/stop the self-rotation
    elif key == glfw.KEY_R:
        View.t_anim = 0                            # restart the pathline animation
        View.t_accum = 0.0
    elif key == glfw.KEY_F:
        View.fade = not View.fade                  # fading tails <-> persistent trails
    elif key in (glfw.KEY_MINUS, glfw.KEY_COMMA):
        View.speed = max(0.25, View.speed * 0.5)   # slow down
    elif key in (glfw.KEY_EQUAL, glfw.KEY_PERIOD):
        View.speed = min(64.0, View.speed * 2.0)   # speed up


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
    path_verts = offsets = lengths = depth_rgb = None
    max_len = 0
    clim = 0.0
    if pl_path and os.path.exists(pl_path):
        print("Loading pathlines:", pl_path)
        path_verts, offsets, lengths, max_len, depths = read_pathlines(
            pl_path, SEED_STRIDE, TIME_STRIDE, EARTH_RADIUS * LINE_LIFT)
        # Per-vertex colour by depth: yellow (surface) -> red (>= clim).
        clim = max(float(np.percentile(depths, DEPTH_CLIM_PCT)), 1.0)
        dnorm = np.clip(depths / clim, 0.0, 1.0)[:, None]
        depth_rgb = ((1.0 - dnorm) * np.array(SURFACE_COLOR, np.float32)
                     + dnorm * np.array(DEEP_COLOR, np.float32)).astype(np.float32)
        print(f"  {len(lengths):,} seeds (every {SEED_STRIDE}), "
              f"{path_verts.shape[0]:,} points, {max_len} animation steps")
        print(f"  depth {depths.min():.0f}-{depths.max():.0f} m; "
              f"colour saturates at p{DEPTH_CLIM_PCT} = {clim:.0f} m")
    else:
        print("No pathlines file — showing globe only. (looked for", pl_path, ")")

    glEnableClientState(GL_VERTEX_ARRAY)

    spin = 0.0
    text_cache = {}
    legend_cache = {}

    # Fade window expressed in TIME: how many subsampled samples span FADE_TIMESTEPS.
    # Alpha ramp indexed by age (index 0 = oldest/transparent, last = newest/opaque).
    fade_pts = max(2, FADE_TIMESTEPS // TIME_STRIDE)
    fade_alpha = np.linspace(0.0, 1.0, fade_pts, dtype=np.float32)
    while not glfw.window_should_close(win):
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

        # --- draw the textured globe (vertex arrays; one glDrawElements) ---
        glEnableClientState(GL_NORMAL_ARRAY)
        glEnableClientState(GL_TEXTURE_COORD_ARRAY)
        glEnable(GL_TEXTURE_2D)
        glVertexPointer(3, GL_FLOAT, 0, verts)
        glNormalPointer(GL_FLOAT, 0, normals)
        glTexCoordPointer(2, GL_FLOAT, 0, uvs)
        glBindTexture(GL_TEXTURE_2D, texture)
        glColor3f(1, 1, 1)
        glDrawElements(GL_TRIANGLES, len(indices), GL_UNSIGNED_INT, indices)

        # --- draw pathlines: growing trails ---
        if path_verts is not None:
            # Only VERTEX_ARRAY should be active, else GL reads normals/uvs
            # (sized for the globe) out of bounds for the larger path buffers.
            glDisableClientState(GL_NORMAL_ARRAY)
            glDisableClientState(GL_TEXTURE_COORD_ARRAY)
            glDisable(GL_TEXTURE_2D)

            t = View.t_anim
            glLineWidth(TRAIL_WIDTH)

            if View.fade:
                # Time-based fade: a point's alpha depends on its age (t - step) in
                # the GLOBAL timeline, so even finished seeds keep fading and vanish
                # once older than the fade_pts window.  RGB comes from depth, so we
                # merge the depth colour with the age alpha into an RGBA slice.
                glEnable(GL_BLEND)
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
                glEnableClientState(GL_COLOR_ARRAY)
                fp1 = fade_pts - 1
                win_start = max(0, t - fp1)                   # oldest step still visible
                for o, L in zip(offsets.tolist(), lengths.tolist()):
                    end = t if t < L else L - 1               # newest existing sample
                    if end < win_start or end - win_start < 1:
                        continue                             # fully faded, or too short
                    idx_end = fp1 - (t - end)                # alpha index for newest point
                    if idx_end < 1:
                        continue                             # newest point already ~invisible
                    idx_start = fp1 - (t - win_start)        # alpha index for oldest point
                    count = end - win_start + 1
                    rgba = np.empty((count, 4), np.float32)
                    rgba[:, :3] = depth_rgb[o + win_start:o + end + 1]     # colour by depth
                    rgba[:, 3] = fade_alpha[idx_start:idx_end + 1]         # alpha by age
                    glVertexPointer(3, GL_FLOAT, 0, path_verts[o + win_start:o + end + 1])
                    glColorPointer(4, GL_FLOAT, 0, rgba)
                    glDrawArrays(GL_LINE_STRIP, 0, count)
                glDisableClientState(GL_COLOR_ARRAY)
                glDisable(GL_BLEND)
            else:
                # Persistent trails: full path 0..t, coloured by depth (yellow->red).
                counts = np.minimum(t + 1, lengths)
                glEnableClientState(GL_COLOR_ARRAY)
                glVertexPointer(3, GL_FLOAT, 0, path_verts)
                glColorPointer(3, GL_FLOAT, 0, depth_rgb)
                for o, c in zip(offsets.tolist(), counts.tolist()):
                    if c >= 2:
                        glDrawArrays(GL_LINE_STRIP, o, c)
                glDisableClientState(GL_COLOR_ARRAY)

            # Current timestep readout, top-left (real MPAS step = t * stride).
            real_t = View.t_anim * TIME_STRIDE
            real_max = (max_len - 1) * TIME_STRIDE
            mode = "fading" if View.fade else "full"
            draw_text_overlay(win,
                              f"timestep  {real_t:,} / {real_max:,}    trails: {mode}"
                              f"    speed {View.speed:g}",
                              text_cache)
            draw_depth_legend(win, clim, legend_cache)

        glfw.swap_buffers(win)
        glfw.poll_events()
        if View.auto_rotate and not View.dragging:
            spin += 0.15          # eastward (real Earth spin): features drift right
        if path_verts is not None and max_len > 1:
            View.t_accum += View.speed          # fractional advance (supports slow speeds)
            if View.t_accum >= max_len:
                View.t_accum = 0.0              # loop the animation
            View.t_anim = int(View.t_accum)

    glfw.terminate()


if __name__ == "__main__":
    main()
