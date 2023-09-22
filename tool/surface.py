#!/usr/bin/python

import sys
import array
import numpy as np
import skimage.measure
import os
import nibabel as nib
import meshio

def readNii(path):
    return nib.load(path).get_fdata()

a = np.array(readNii(sys.argv[1]))
nx, ny, nz = np.shape(a)
node = np.zeros_like(a)
cnt = 0
w = 0, 1
for dx in w:
    for dy in w:
        for dz in w:
            q = 1 / ((1 + abs(dx)) * (1 + abs(dy)) * (1 + abs(dz)))
            print(dx, dy, dz, q)
            node = node + q * np.roll(a, (dx, dy, dz))
            cnt += q
node /= cnt
print(np.std(node), np.std(a), np.shape(node), np.shape(a))
verts, faces, normals, values = skimage.measure.marching_cubes(node,
                                                               0.9,
                                                               spacing=(1/nx, 1/ny, 1/nz))
mesh = meshio.Mesh(
    verts,
    (("triangle", faces),)
)
mesh.write("surface.vtk")

path = "surface"
xyz_path = "%s.xyz.raw" % path
topo_path = "%s.topo.raw" % path
xdmf_path = "%s.xdmf2" % path
verts.tofile(xyz_path)
faces.tofile(topo_path)
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
