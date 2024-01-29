import SimpleITK as sitk
import vtk
from vtk.util import numpy_support


def sitk_to_vtk(sitk_image):
    # Convert SimpleITK image to VTK image for processing
    sitk_array = sitk.GetArrayFromImage(sitk_image)
    vtk_data_array = numpy_support.numpy_to_vtk(num_array=sitk_array.ravel(), deep=True, array_type=vtk.VTK_FLOAT)

    vtk_image = vtk.vtkImageData()
    vtk_image.SetDimensions(sitk_image.GetSize())
    vtk_image.GetPointData().SetScalars(vtk_data_array)
    return vtk_image


# Read the NIfTI file
nii_file = 'output/9E27833A/203_bone.nii.gz'
sitk_image = sitk.ReadImage(nii_file)

# Apply thresholding to isolate the skull
lower_threshold = 1000  # Adjust these values based on your data
upper_threshold = 2000
bone_image = sitk.BinaryThreshold(sitk_image, lower_threshold, upper_threshold, 1, 0)

# Convert to VTK image
vtk_image = sitk_to_vtk(bone_image)

# Create a 3D surface (mesh) using VTK
contour_filter = vtk.vtkContourFilter()
contour_filter.SetInputData(vtk_image)
contour_filter.SetValue(0, 0.5)  # Threshold for contouring

# Optional: Smooth the mesh
smoother = vtk.vtkWindowedSincPolyDataFilter()
smoother.SetInputConnection(contour_filter.GetOutputPort())
smoother.SetNumberOfIterations(15)
smoother.NonManifoldSmoothingOn()
smoother.NormalizeCoordinatesOn()
smoother.Update()

# Export to STL or PLY
# stl_writer = vtk.vtkSTLWriter()
# stl_writer.SetFileName('skull_model.stl')
# stl_writer.SetInputConnection(smoother.GetOutputPort())
# stl_writer.Write()

ply_writer = vtk.vtkPLYWriter()
ply_writer.SetFileName('skull_model.ply')
ply_writer.SetInputConnection(smoother.GetOutputPort())
ply_writer.Write()