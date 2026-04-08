"""
Read and reassemble streamlines from trace_outputs/traces_<gid>.bin files.

Input binary format per segment:
  [int32 pid][int32 gid][int32 sid][int32 nPts]
  [float64 x0][float64 y0][float64 z0] ...

Output binary format (streamlines.bin):
  [int32 nStreamlines]
  For each streamline (sorted by pid):
    [int32 pid][int32 nPts][float64 x0][float64 y0][float64 z0] ...

Usage:
  python read_traces.py [nblocks] [trace_dir]
  python read_traces.py 16 trace_outputs
"""

import struct
import sys
import os
from collections import defaultdict


def read_trace_file(filename):
    """Read one traces_<gid>.bin and return list of segment dicts."""
    segments = []
    with open(filename, 'rb') as f:
        while True:
            header = f.read(16)          # 4 × int32
            if len(header) < 16:
                break
            pid, gid, sid, nPts = struct.unpack('iiii', header)
            raw = f.read(nPts * 24)      # nPts × 3 × float64
            if len(raw) < nPts * 24:
                print(f'Warning: truncated data in {filename} (pid={pid} sid={sid})')
                break
            coords = list(struct.iter_unpack('ddd', raw))
            segments.append({'pid': pid, 'gid': gid, 'sid': sid, 'coords': coords})
    return segments


def reassemble_streamlines(all_segments):
    """
    Group segments by pid, sort by sid, concatenate into full streamlines.
    The boundary point shared between consecutive segments is deduplicated
    (segment N's last point == segment N+1's first point).
    Returns: dict  pid -> list of (x, y, z)
    """
    by_pid = defaultdict(list)
    for seg in all_segments:
        by_pid[seg['pid']].append(seg)

    streamlines = {}
    for pid, segs in sorted(by_pid.items()):
        segs.sort(key=lambda s: s['sid'])
        pts = []
        for i, seg in enumerate(segs):
            # Skip first point of each segment after the first —
            # it duplicates the last point of the previous segment.
            start = 0 if i == 0 else 1
            pts.extend(seg['coords'][start:])
        streamlines[pid] = pts
    return streamlines


def print_summary(streamlines):
    total_pts = sum(len(pts) for pts in streamlines.values())
    print(f'Streamlines : {len(streamlines)}')
    print(f'Total points: {total_pts}')
    print()
    for pid, pts in sorted(streamlines.items()):
        print(f'  pid={pid:4d}  {len(pts):6d} pts   '
              f'start_radius={(pts[0][0]**2 + pts[0][1]**2 + pts[0][2]**2)**0.5:.3f}  '
              f'end_radius={(pts[-1][0]**2 + pts[-1][1]**2 + pts[-1][2]**2)**0.5:.3f}  '
              f'start=({pts[0][0]:.3f}, {pts[0][1]:.3f}, {pts[0][2]:.3f})  '
              f'end=({pts[-1][0]:.3f}, {pts[-1][1]:.3f}, {pts[-1][2]:.3f})')


def save_as_binary(streamlines, output_path):
    """
    Save reassembled streamlines to a binary file.

    Format:
      [int32 nStreamlines]
      For each streamline (sorted by pid):
        [int32 pid][int32 nPts]
        [float64 x0][float64 y0][float64 z0]
        ...
    """
    sorted_pids = sorted(streamlines.keys())
    with open(output_path, 'wb') as f:
        f.write(struct.pack('i', len(sorted_pids)))
        for pid in sorted_pids:
            pts = streamlines[pid]
            f.write(struct.pack('ii', pid, len(pts)))
            for x, y, z in pts:
                f.write(struct.pack('ddd', x, y, z))
    print(f'Saved {len(sorted_pids)} streamlines to {output_path}')


def main():
    nblocks   = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    trace_dir = sys.argv[2]      if len(sys.argv) > 2 else 'streamlines'

    all_segments = []
    for gid in range(nblocks):
        path = os.path.join(trace_dir, f'{gid}.bin')
        if not os.path.exists(path):
            print(f'Warning: {path} not found, skipping')
            continue
        segs = read_trace_file(path)
        print(f'gid {gid:2d}: {len(segs)} segments')
        all_segments.extend(segs)

    print(f'\nTotal segments across all blocks: {len(all_segments)}\n')

    streamlines = reassemble_streamlines(all_segments)
    print_summary(streamlines)

    out_bin = os.path.join(trace_dir, trace_dir+'.bin')
    save_as_binary(streamlines, out_bin)


if __name__ == '__main__':
    main()