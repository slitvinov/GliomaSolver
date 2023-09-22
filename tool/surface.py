#!/usr/bin/python

import sys
import array
import numpy as np
import skimage.measure
import os
import nibabel as nib

def readNii(path):
    return nib.load(path).get_fdata()

a = readNii(sys.argv[1])
sx, sy, sz = 1, 1, 1
nx, ny, nz = np.shape(a)
dx = 1.0 / nx / sx
dy = 1.0 / ny / sy
dz = 1.0 / nz / sz
verts, faces, normals, values = skimage.measure.marching_cubes(a[::sx, ::sy, ::sz],
                                                               0.5,
                                                               spacing=(dx, dy,
                                                                        dz))
path = "surface"
xyz_path = "%s.xyz.raw" % path
topo_path = "%s.topo.raw" % path
xdmf_path = "%s.xdmf2" % path
xyz = np.memmap(xyz_path, np.dtype("<f4"), "w+", shape=(len(verts), 3), order='C')
topo = np.memmap(topo_path, np.dtype("<i4"), "w+", shape=(len(faces), 3), order='C')
np.copyto(xyz, verts)
np.copyto(topo, faces, 'no')
with open(xdmf_path, "w") as f:
    f.write("""\
<Xdmf
    Version="2">
  <Domain>
    <Grid>
      <Topology
	  TopologyType="Triangle"
	  Dimensions="%d">
	<DataItem
	    Dimensions="%d 3"
	    NumberType="Int"
	    Format="Binary">
          %s
	</DataItem>
      </Topology>
      <Geometry>
	<DataItem
	    Dimensions="%d 3"
	    Format="Binary">
          %s
	</DataItem>
      </Geometry>
    </Grid>
  </Domain>
</Xdmf>
""" % (len(faces), len(faces), os.path.basename(topo_path), len(verts),
       os.path.basename(xyz_path)))
