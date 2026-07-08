"""
Generate a METIS partition file from an MPAS-O NetCDF mesh.

Usage:
    python3 gen_partition.py <mesh.nc> <n_parts> <output_prefix>

Examples:
    # Lowres, 256 blocks (CPU jobs)
    python3 gen_partition.py /global/cfs/cdirs/m4259/yqiu/MOPS_ab/climatology/ocean.EC30to60E2r2.210210.nc 256 lowres_20y

    # Lowres, 16 blocks (GPU jobs, 1 task per node)
    python3 gen_partition.py /global/cfs/cdirs/m4259/yqiu/MOPS_ab/climatology/ocean.EC30to60E2r2.210210.nc 16 lowres_20y

    # Highres, 256 blocks
    python3 gen_partition.py 20210421_sim7_CORE_18to6v3.mpaso.rst.0015-01-01_00000.nc 256 highres_6y

    # Highres, 16 blocks
    python3 gen_partition.py 20210421_sim7_CORE_18to6v3.mpaso.rst.0015-01-01_00000.nc 16 highres_6y

Output files:
    <prefix>.info           METIS graph file (adjacency list, 1-indexed)
    <prefix>.info.part.<N>  Partition assignment: one integer per line (0..N-1),
                            one line per cell. This is what ParaFlow reads via
                            graph_partition_indices in the YAML config.
"""

import os
import sys
import netCDF4 as nc
import pymetis


def build_adjacency(mesh_file):
    """
    Read cell-to-cell adjacency from MPAS mesh NetCDF.
    Returns (adjacency, n_cells) where adjacency[i] is a list of 0-indexed
    neighbor cell IDs for cell i.
    """
    print(f"Reading mesh: {mesh_file}")
    ds = nc.Dataset(mesh_file, 'r')
    cells_on_cell   = ds.variables['cellsOnCell'][:]   # (nCells, maxEdges), 1-indexed; 0 = no neighbor
    n_edges_on_cell = ds.variables['nEdgesOnCell'][:]  # (nCells,)
    ds.close()

    n_cells = cells_on_cell.shape[0]
    print(f"  {n_cells} cells, max {cells_on_cell.shape[1]} edges per cell")

    adjacency = []
    for i in range(n_cells):
        n_edges = int(n_edges_on_cell[i])
        neighbors = [int(cells_on_cell[i, j]) - 1   # convert to 0-indexed
                     for j in range(n_edges)
                     if cells_on_cell[i, j] > 0]    # 0 means boundary padding
        adjacency.append(neighbors)

    return adjacency, n_cells


def write_metis_graph(adjacency, n_cells, graph_file):
    """Write METIS graph.info format for reference / gpmetis compatibility."""
    n_edges = sum(len(adj) for adj in adjacency) // 2
    with open(graph_file, 'w') as f:
        f.write(f"{n_cells} {n_edges}\n")
        for neighbors in adjacency:
            f.write(" ".join(str(n + 1) for n in neighbors) + "\n")
    print(f"  Wrote METIS graph file: {graph_file}  ({n_cells} cells, {n_edges} edges)")


def partition_and_write(adjacency, n_parts, part_file):
    """Partition with pymetis and write one-integer-per-line output."""
    print(f"  Partitioning into {n_parts} parts...")
    _, membership = pymetis.part_graph(n_parts, adjacency=adjacency)

    counts = [0] * n_parts
    for p in membership:
        counts[p] += 1
    print(f"  Cells per partition: min={min(counts)}  max={max(counts)}  avg={len(membership)//n_parts}")

    with open(part_file, 'w') as f:
        for part_id in membership:
            f.write(f"{part_id}\n")
    print(f"  Wrote partition file: {part_file}")


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    mesh_file = sys.argv[1]
    n_parts   = int(sys.argv[2])
    prefix    = sys.argv[3]

    out_dir    = os.path.join(os.path.dirname(os.path.abspath(__file__)), "partition")
    os.makedirs(out_dir, exist_ok=True)
    graph_file = os.path.join(out_dir, f"{prefix}.info")
    part_file  = os.path.join(out_dir, f"{prefix}.info.part.{n_parts}")

    adjacency, n_cells = build_adjacency(mesh_file)
    write_metis_graph(adjacency, n_cells, graph_file)
    partition_and_write(adjacency, n_parts, part_file)
    print("Done.")


if __name__ == "__main__":
    main()
