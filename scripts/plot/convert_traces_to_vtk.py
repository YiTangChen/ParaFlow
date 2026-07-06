#!/usr/bin/env python3
"""
Convert MPAS-O particle trace binary files to VTK PolyData (.vtp)
for visualization in ParaView.

Reads the streamlines.bin file produced by read_traces.py.

streamlines.bin format:
  int32 nStreamlines
  For each streamline:
    int32 pid
    int32 nPts
    nPts * [float64 x, float64 y, float64 z]

Usage:
    python3 convert_traces_to_vtk.py [streamlines.bin] [output.vtp]

Defaults:
    streamlines.bin : trace_outputs/streamlines.bin
    output          : traces.vtp
"""

import sys
import os
import struct
import numpy as np


def read_streamlines_bin(filepath):
    """
    Read streamlines.bin produced by read_traces.py.
    Returns list of (pid, Nx3 numpy array) tuples.
    """
    streamlines = []
    with open(filepath, 'rb') as f:
        data = f.read()
    offset = 0
    total_size = len(data)

    if total_size < 4:
        raise ValueError(f"File too small: {filepath}")
    n_streamlines = struct.unpack_from('<i', data, offset)[0]
    offset += 4
    print(f"  nStreamlines = {n_streamlines}")

    for _ in range(n_streamlines):
        if offset + 8 > total_size:
            print("  Warning: unexpected end of file")
            break
        pid, nPts = struct.unpack_from('<ii', data, offset)
        offset += 8
        byte_count = nPts * 3 * 8
        if offset + byte_count > total_size:
            print(f"  Warning: truncated data for pid={pid}, stopping")
            break
        coords = np.frombuffer(data, dtype='<f8', count=nPts * 3, offset=offset)
        coords = coords.reshape(nPts, 3)
        streamlines.append((pid, coords))
        offset += byte_count

    return streamlines


def read_all_traces(trace_dir):
    """Read streamlines.bin from trace_dir, return list of (pid, Nx3) tuples."""
    filepath = os.path.join(trace_dir, 'streamlines.bin')
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"streamlines.bin not found: {filepath}\n"
                                f"Run read_traces.py first to generate it.")
    print(f"  Reading {filepath} ...")
    return read_streamlines_bin(filepath)


def write_vtp(streamlines, output_path):
    """Write streamlines as VTK XML PolyData (.vtp) with polylines."""
    try:
        import vtk
        _write_vtp_vtk(streamlines, output_path)
    except ImportError:
        print("  vtk module not found, writing raw XML VTK format...")
        _write_vtp_xml(streamlines, output_path)


def _write_vtp_vtk(streamlines, output_path):
    """Write using vtk Python library."""
    import vtk

    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    pid_arr = vtk.vtkIntArray()
    pid_arr.SetName("ParticleID")
    pid_arr.SetNumberOfComponents(1)

    pt_offset = 0
    for pid, coords in streamlines:
        nPts = len(coords)
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(nPts)
        for i, (x, y, z) in enumerate(coords):
            points.InsertNextPoint(x, y, z)
            polyline.GetPointIds().SetId(i, pt_offset + i)
        lines.InsertNextCell(polyline)
        pid_arr.InsertNextValue(pid)
        pt_offset += nPts

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetCellData().AddArray(pid_arr)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_path)
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    writer.Write()
    print(f"  Wrote {output_path} using vtk library")


def _write_vtp_xml(streamlines, output_path):
    """Fallback: write VTK XML PolyData manually (binary base64 encoded)."""
    import base64

    coord_arrays = [coords for _, coords in streamlines]
    all_pts = np.concatenate(coord_arrays, axis=0).astype(np.float32)
    nTotalPts = len(all_pts)

    offsets = []
    connectivity = []
    pid_per_cell = []
    cumulative = 0
    for pid, coords in streamlines:
        n = len(coords)
        connectivity.extend(range(cumulative, cumulative + n))
        cumulative += n
        offsets.append(cumulative)
        pid_per_cell.append(pid)

    conn_arr = np.array(connectivity, dtype=np.int32)
    offs_arr = np.array(offsets, dtype=np.int32)
    pid_arr  = np.array(pid_per_cell, dtype=np.int32)

    def encode(arr):
        raw = arr.tobytes()
        header = struct.pack('<Q', len(raw))
        return base64.b64encode(header + raw).decode('ascii')

    pts_flat = all_pts.flatten()

    xml_lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt64">',
        '  <PolyData>',
        f'    <Piece NumberOfPoints="{nTotalPts}" NumberOfVerts="0" '
        f'NumberOfLines="{len(streamlines)}" NumberOfStrips="0" NumberOfPolys="0">',
        '      <PointData/>',
        '      <CellData>',
        f'        <DataArray type="Int32" Name="ParticleID" format="binary">',
        f'          {encode(pid_arr)}',
        '        </DataArray>',
        '      </CellData>',
        '      <Points>',
        '        <DataArray type="Float32" NumberOfComponents="3" format="binary">',
        f'          {encode(pts_flat)}',
        '        </DataArray>',
        '      </Points>',
        '      <Lines>',
        '        <DataArray type="Int32" Name="connectivity" format="binary">',
        f'          {encode(conn_arr)}',
        '        </DataArray>',
        '        <DataArray type="Int32" Name="offsets" format="binary">',
        f'          {encode(offs_arr)}',
        '        </DataArray>',
        '      </Lines>',
        '    </Piece>',
        '  </PolyData>',
        '</VTKFile>',
    ]
    with open(output_path, 'w') as f:
        f.write('\n'.join(xml_lines))
    print(f"  Wrote {output_path} using fallback XML writer")


def write_netcdf(streamlines, output_path):
    """
    Write streamlines as NetCDF4 (CF-conventions style ragged array).
    VTK format is recommended for direct line display in ParaView.
    """
    try:
        import netCDF4 as nc
    except ImportError:
        print("netCDF4 module not found, skipping NetCDF output")
        return

    coord_arrays = [coords for _, coords in streamlines]
    all_pts = np.concatenate(coord_arrays, axis=0)
    lengths  = np.array([len(c) for c in coord_arrays], dtype=np.int32)
    offsets  = np.concatenate([[0], np.cumsum(lengths[:-1])]).astype(np.int32)
    pids     = np.array([pid for pid, _ in streamlines], dtype=np.int32)
    pid_per_pt = np.repeat(pids, lengths)

    ds = nc.Dataset(output_path, 'w')
    ds.createDimension('nPoints', len(all_pts))
    ds.createDimension('nTraces', len(streamlines))

    vx = ds.createVariable('x', 'f8', ('nPoints',))
    vy = ds.createVariable('y', 'f8', ('nPoints',))
    vz = ds.createVariable('z', 'f8', ('nPoints',))
    vx[:] = all_pts[:, 0];  vx.units = 'meters'
    vy[:] = all_pts[:, 1];  vy.units = 'meters'
    vz[:] = all_pts[:, 2];  vz.units = 'meters'

    vo  = ds.createVariable('trace_offset', 'i4', ('nTraces',))
    vl  = ds.createVariable('trace_length', 'i4', ('nTraces',))
    vp  = ds.createVariable('pid',          'i4', ('nTraces',))
    vpt = ds.createVariable('pid_per_pt',   'i4', ('nPoints',))
    vo[:]  = offsets
    vl[:]  = lengths
    vp[:]  = pids
    vpt[:] = pid_per_pt

    ds.description = 'MPAS-O particle traces (Cartesian XYZ in meters)'
    ds.close()
    print(f"  Wrote {output_path}")


def main():
    # Accept either a streamlines.bin path or a trace_dir
    if len(sys.argv) > 1 and sys.argv[1].endswith('.bin'):
        bin_path = sys.argv[1]
        trace_dir = os.path.dirname(bin_path) or '.'
    else:
        trace_dir = sys.argv[1] if len(sys.argv) > 1 else 'trace_outputs'
        bin_path  = None   # read_all_traces will find streamlines.bin in trace_dir

    output = sys.argv[2] if len(sys.argv) > 2 else 'traces.vtp'

    if bin_path:
        print(f"Reading: {bin_path}")
        streamlines = read_streamlines_bin(bin_path)
    else:
        print(f"Reading from: {trace_dir}")
        streamlines = read_all_traces(trace_dir)

    total_pts = sum(len(c) for _, c in streamlines)
    print(f"Total: {len(streamlines)} streamlines, {total_pts:,} points")

    ext = os.path.splitext(output)[1].lower()
    if ext in ('.nc', '.netcdf'):
        print(f"Writing NetCDF: {output}")
        write_netcdf(streamlines, output)
    else:
        print(f"Writing VTK PolyData: {output}")
        write_vtp(streamlines, output)

    print("Done.")


if __name__ == '__main__':
    main()
