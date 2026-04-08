#!/usr/bin/env python3
"""
Parse memory output from ParaFlow MPI runs.

Expected log line format (written to stderr by each MPI rank):
  MEM_AFTER_LOAD rank=0 gid=3 VmRSS:  1234567 kB
  MEM_PEAK       rank=0 gid=3 VmHWM:  1234567 kB

Usage:
  # single file
  python3 parse_mem.py run.log

  # compare before / after
  python3 parse_mem.py before.log after.log
"""

import re
import sys
from pathlib import Path

# matches both MEM_AFTER_LOAD and MEM_PEAK lines
PATTERN = re.compile(
    r"(MEM_\w+)\s+rank=(\d+)\s+gid=(\d+)\s+\w+:\s+(\d+)\s+kB"
)


def parse_log(path: str) -> dict:
    """Return {tag: {gid: kb}} from a log file."""
    result = {}  # type: dict
    with open(path) as f:
        for line in f:
            m = PATTERN.search(line)
            if not m:
                continue
            tag, rank, gid, kb = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
            result.setdefault(tag, {})[gid] = kb
    return result


def summarize(label: str, data: dict):
    if not data:
        print(f"  [{label}] No MEM lines found.")
        return
    for tag, gid_map in sorted(data.items()):
        values = list(gid_map.values())
        n      = len(values)
        avg    = sum(values) / n
        mn     = min(values)
        mx     = max(values)
        total  = sum(values)
        print(f"\n  [{label}] {tag}  ({n} blocks)")
        print(f"    avg   : {avg/1024:>10.1f} MB = {avg/1024/1024:>10.1f} GB")
        print(f"    min   : {mn /1024:>10.1f} MB = {mn /1024/1024:>10.1f} GB  (gid {min(gid_map, key=gid_map.get)})")
        print(f"    max   : {mx /1024:>10.1f} MB = {mx /1024/1024:>10.1f} GB  (gid {max(gid_map, key=gid_map.get)})")
        print(f"    total : {total/1024:>10.1f} MB = {total/1024/1024:>10.1f} GB  (sum across all blocks)")


def compare(before_data: dict, after_data: dict):
    all_tags = set(before_data) | set(after_data)
    for tag in sorted(all_tags):
        b = before_data.get(tag, {})
        a = after_data.get(tag, {})
        common_gids = set(b) & set(a)
        if not common_gids:
            print(f"\n  [{tag}] No common gids to compare.")
            continue

        b_vals = [b[g] for g in common_gids]
        a_vals = [a[g] for g in common_gids]
        b_avg  = sum(b_vals) / len(b_vals)
        a_avg  = sum(a_vals) / len(a_vals)
        b_max  = max(b_vals)
        a_max  = max(a_vals)
        diff_avg = b_avg - a_avg
        diff_max = b_max - a_max
        pct_avg  = diff_avg / b_avg * 100 if b_avg else 0
        pct_max  = diff_max / b_max * 100 if b_max else 0

        print(f"\n  {tag}  ({len(common_gids)} common blocks)")
        print(f"  {'':20s} {'before':>12} {'after':>12} {'saved':>12} {'reduction':>10}")
        print(f"  {'avg per block':20s} {b_avg/1024:>10.1f}MB {a_avg/1024:>10.1f}MB "
              f"{diff_avg/1024:>10.1f}MB {pct_avg:>9.1f}%")
        print(f"  {'max block':20s} {b_max/1024:>10.1f}MB {a_max/1024:>10.1f}MB "
              f"{diff_max/1024:>10.1f}MB {pct_max:>9.1f}%")

        # per-gid breakdown
        print(f"\n  {'gid':>6} {'before(MB)':>12} {'after(MB)':>12} {'diff(MB)':>10} {'reduction':>10}")
        for g in sorted(common_gids):
            bv, av = b[g], a[g]
            d  = bv - av
            pct = d / bv * 100 if bv else 0
            print(f"  {g:>6} {bv/1024:>12.1f} {av/1024:>12.1f} {d/1024:>10.1f} {pct:>9.1f}%")


def main():
    if len(sys.argv) == 2:
        path = sys.argv[1]
        print(f"=== {Path(path).name} ===")
        data = parse_log(path)
        summarize(Path(path).name, data)

    elif len(sys.argv) == 3:
        before_path, after_path = sys.argv[1], sys.argv[2]
        print(f"=== Comparing: {Path(before_path).name}  vs  {Path(after_path).name} ===")
        before_data = parse_log(before_path)
        after_data  = parse_log(after_path)

        print("\n--- BEFORE ---")
        summarize("before", before_data)
        print("\n--- AFTER ---")
        summarize("after", after_data)
        print("\n--- DIFF ---")
        compare(before_data, after_data)

    else:
        print("Usage:")
        print("  python3 parse_mem.py <run.log>")
        print("  python3 parse_mem.py <before.log> <after.log>")
        sys.exit(1)


if __name__ == "__main__":
    main()
