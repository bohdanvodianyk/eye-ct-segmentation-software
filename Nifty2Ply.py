import nibabel as nib
import numpy as np
from stl import mesh
from skimage import measure

file_path = 'test-model/labels/final/5_bone_05.nii.gz'

nifti_file = nib.load(file_path)
np_array = nifti_file.get_fdata()

verts, faces, normals, values = measure.marching_cubes(np_array, 0)

obj_3d = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))

for i, f in enumerate(faces):
    obj_3d.vectors[i] = verts[f]

obj_3d.save('segmentation.stl')
