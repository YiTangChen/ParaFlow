#!/usr/bin/env python3
"""
Parse memory output from ParaFlow MPI runs.

Supported log line formats (written to stderr by each MPI rank):
  MEM_ANALYTICAL rank=0 gid=3 grid_bytes=16965796 solution_bytes=47600640 total_bytes=64566436
  MEM_DELTA      rank=0 gid=3 before_kb=95892 after_kb=392584 delta_kb=296692
  MEM_PEAK       rank=0 gid=3 vmhwm_kb=2620692

Legacy format:
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

# Matches the older "VmRSS: ... kB" / "VmHWM: ... kB" format.
LEGACY_PATTERN = re.compile(
    r"(MEM_\w+)\s+rank=(\d+)\s+gid=(\d+)\s+(\w+):\s+(\d+)\s+kB"
)

# Matches newer MEM_* key=value lines.
MEM_LINE_PATTERN = re.compile(r"\b(MEM_\w+)\b\s+(.*)")
KV_PATTERN = re.compile(r"(\w+)=([+-]?\d+(?:\.\d+)?)")


def parse_number(value: str):
    if "." in value:
        return float(value)
    return int(value)


def parse_log(path: str) -> dict:
    """Return {tag: {metric: {gid: value}}} from a log file."""
    result = {}  # type: dict
    with open(path) as f:
        for line in f:
            legacy = LEGACY_PATTERN.search(line)
            if legacy:
                tag = legacy.group(1)
                gid = int(legacy.group(3))
                metric = f"{legacy.group(4).lower()}_kb"
                kb = int(legacy.group(5))
                result.setdefault(tag, {}).setdefault(metric, {})[gid] = kb
                continue

            mem_line = MEM_LINE_PATTERN.search(line)
            if not mem_line:
                continue

            tag = mem_line.group(1)
            fields = dict(KV_PATTERN.findall(mem_line.group(2)))
            if "gid" not in fields:
                continue

            gid = int(fields["gid"])
            for metric, value in fields.items():
                if metric in ("rank", "gid"):
                    continue
                result.setdefault(tag, {}).setdefault(metric, {})[gid] = parse_number(value)
    return result


def format_value(value, metric: str) -> str:
    if metric.endswith("_bytes"):
        return f"{value/1024/1024:>10.1f} MB = {value/1024/1024/1024:>10.1f} GB"
    if metric.endswith("_kb"):
        return f"{value/1024:>10.1f} MB = {value/1024/1024:>10.1f} GB"
    return f"{value:>10.1f}"


def format_delta(value, metric: str) -> str:
    if metric.endswith("_bytes"):
        return f"{value/1024/1024:>10.1f} MB"
    if metric.endswith("_kb"):
        return f"{value/1024:>10.1f} MB"
    return f"{value:>10.1f}"


def summarize(label: str, data: dict):
    if not data:
        print(f"  [{label}] No MEM lines found.")
        return
    for tag, metrics in sorted(data.items()):
        print(f"\n  [{label}] {tag}")
        for metric, gid_map in sorted(metrics.items()):
            values = list(gid_map.values())
            n      = len(values)
            avg    = sum(values) / n
            mn     = min(values)
            mx     = max(values)
            total  = sum(values)
            print(f"    {metric}  ({n} blocks)")
            print(f"      avg   : {format_value(avg, metric)}")
            print(f"      min   : {format_value(mn, metric)}  (gid {min(gid_map, key=gid_map.get)})")
            print(f"      max   : {format_value(mx, metric)}  (gid {max(gid_map, key=gid_map.get)})")
            print(f"      total : {format_value(total, metric)}  (sum across all blocks)")


def compare(before_data: dict, after_data: dict):
    all_tags = set(before_data) | set(after_data)
    for tag in sorted(all_tags):
        before_metrics = before_data.get(tag, {})
        after_metrics = after_data.get(tag, {})
        all_metrics = set(before_metrics) | set(after_metrics)
        if not all_metrics:
            print(f"\n  [{tag}] No metrics to compare.")
            continue

        for metric in sorted(all_metrics):
            b = before_metrics.get(metric, {})
            a = after_metrics.get(metric, {})
            common_gids = set(b) & set(a)
            if not common_gids:
                print(f"\n  [{tag} {metric}] No common gids to compare.")
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

            print(f"\n  {tag} {metric}  ({len(common_gids)} common blocks)")
            print(f"  {'':20s} {'before':>14} {'after':>14} {'saved':>14} {'reduction':>10}")
            print(f"  {'avg per block':20s} {format_delta(b_avg, metric)} {format_delta(a_avg, metric)} "
                  f"{format_delta(diff_avg, metric)} {pct_avg:>9.1f}%")
            print(f"  {'max block':20s} {format_delta(b_max, metric)} {format_delta(a_max, metric)} "
                  f"{format_delta(diff_max, metric)} {pct_max:>9.1f}%")

            # per-gid breakdown
            print(f"\n  {'gid':>6} {'before(MB)':>12} {'after(MB)':>12} {'diff(MB)':>10} {'reduction':>10}")
            for g in sorted(common_gids):
                bv, av = b[g], a[g]
                d  = bv - av
                pct = d / bv * 100 if bv else 0
                print(f"  {g:>6} {format_delta(bv, metric)} {format_delta(av, metric)} "
                      f"{format_delta(d, metric)} {pct:>9.1f}%")


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
