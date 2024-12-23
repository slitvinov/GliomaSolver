#!/usr/bin/python

import sys
import array
import numpy as np
import skimage.measure
import scipy.ndimage
import os
import nibabel as nib


def readNii(path):
    return nib.load(path).get_fdata()


path = sys.argv[1]
try:
    a = readNii(path)
except FileNotFoundError as e:
    sys.stderr.write("surface.py: error: %s\n" % e)
    sys.exit(1)
except nib.filebasedimages.ImageFileError as e:
    sys.stderr.write("surface.py: error: %s\n" % e)
    sys.exit(1)

sx, sy, sz = 1, 1, 1
nx, ny, nz = np.shape(a)
node = np.zeros_like(a)
scipy.ndimage.gaussian_filter(a, sigma=1, output=node)
verts, faces, normals, values = skimage.measure.marching_cubes(
    node[::sx, ::sy, ::sz], 0.5, spacing=(sx / nx, sy / ny, sz / nz))
path = "surface"
xyz_path = "%s.xyz.raw" % path
topo_path = "%s.topo.raw" % path
xdmf_path = "%s.xdmf2" % path
verts.tofile(xyz_path)
faces.tofile(topo_path)
sys.stderr.write("surface.py: %s\n" % xdmf_path)
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
	    Precision="%d"
	    Format="Binary">
          %s
	</DataItem>
      </Geometry>
    </Grid>
  </Domain>
</Xdmf>
""" % (len(faces), len(faces), os.path.basename(topo_path), len(verts),
       verts.itemsize, os.path.basename(xyz_path)))
