#!/bin/bash

# python3 read_traces.py 256 streamlines
# python3 plot_traces.py streamlines/streamlines.bin drawSubdomain streamlines_map.png

python3 read_traces.py 256 pathlines
python3 plot_traces.py pathlines/pathlines.bin drawSubdomain pathlines_map.png

