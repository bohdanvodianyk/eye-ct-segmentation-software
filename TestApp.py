import tkinter as tk
from tkinter import filedialog
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Slider
import SimpleITK as sitk
import vtk
from vtk.util import numpy_support
from tkinter import font as tkFont

sitk_image = None
lower_threshold = None
upper_threshold = None


def apply_threshold_and_normalize(slice, lower_threshold, upper_threshold):
    thresholded_slice = np.copy(slice)
    thresholded_slice[slice < lower_threshold] = 0
    thresholded_slice[slice > upper_threshold] = 0
    valid_mask = (thresholded_slice > 0)
    if np.any(valid_mask):
        normalized_slice = np.zeros_like(thresholded_slice, dtype=np.float32)
        normalized_slice[valid_mask] = (thresholded_slice[valid_mask] - np.min(slice)) / (np.max(slice) - np.min(slice))
        normalized_slice *= 255
        return normalized_slice
    else:
        return thresholded_slice


def load_and_visualize_nifti(visualization_window, canvas, ax, fig):
    global sitk_image
    file_path = filedialog.askopenfilename(filetypes=[("NIfTI files", "*.nii;*.nii.gz")])
    if not file_path:
        return

    nifti_image = nib.load(file_path)
    sitk_image = sitk.ReadImage(file_path)
    image_array = nifti_image.get_fdata()

    slice_idx = 0
    lower_threshold, upper_threshold = 0, np.max(image_array)
    displayed_image = apply_threshold_and_normalize(image_array[:, :, slice_idx], lower_threshold, upper_threshold)
    img = ax.imshow(displayed_image, cmap="gray", vmin=0, vmax=255)

    ax_slice = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slice_slider = Slider(ax_slice, 'Slice', 0, image_array.shape[2] - 1, valinit=0, valstep=1)

    ax_lower = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    lower_slider = Slider(ax_lower, 'Lower Threshold', np.min(image_array), np.max(image_array), valinit=0)

    ax_upper = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    upper_slider = Slider(ax_upper, 'Upper Threshold', np.min(image_array), np.max(image_array),
                          valinit=upper_threshold)

    def update(val):
        global lower_threshold, upper_threshold
        slice_idx = int(slice_slider.val)
        lower_threshold = lower_slider.val
        upper_threshold = upper_slider.val
        normalized_image = apply_threshold_and_normalize(image_array[:, :, slice_idx], lower_threshold, upper_threshold)
        img.set_data(normalized_image)
        canvas.draw_idle()

    slice_slider.on_changed(update)
    lower_slider.on_changed(update)
    upper_slider.on_changed(update)

    canvas.draw()


def sitk_to_vtk(sitk_image):
    # Convert SimpleITK image to VTK image for processing
    sitk_array = sitk.GetArrayFromImage(sitk_image)
    vtk_data_array = numpy_support.numpy_to_vtk(num_array=sitk_array.ravel(), deep=True, array_type=vtk.VTK_FLOAT)

    vtk_image = vtk.vtkImageData()
    vtk_image.SetDimensions(sitk_image.GetSize())
    vtk_image.GetPointData().SetScalars(vtk_data_array)
    return vtk_image


def perform_3d_visualization():
    global sitk_image, lower_threshold, upper_threshold
    print(lower_threshold, upper_threshold)
    bone_image = sitk.BinaryThreshold(sitk_image, lower_threshold, upper_threshold, 1, 0)
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

    ply_writer = vtk.vtkPLYWriter()
    ply_writer.SetFileName('test.ply')
    ply_writer.SetInputConnection(smoother.GetOutputPort())
    ply_writer.Write()
    # Implement 3D visualization using VTK here
    print("3D visualization to be implemented")


# Create the main visualization window
visualization_window = tk.Tk()
visualization_window.title("3D Visualization")
helv36 = tkFont.Font(family='Helvetica', size=20, weight='bold')

fig, ax = plt.subplots()
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.3, top=0.9)  # Adjust subplot size

canvas = FigureCanvasTkAgg(fig, master=visualization_window)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Create frames for button rows
button_frame1 = tk.Frame(visualization_window)
button_frame2 = tk.Frame(visualization_window)
button_frame3 = tk.Frame(visualization_window)

# Create buttons
load_button = tk.Button(button_frame1, text="Open NIfTI File", command=lambda: load_and_visualize_nifti(visualization_window, canvas, ax, fig))
load_button['font'] = helv36
dicom_load = tk.Button(button_frame1, text="Open DICOM folder")
dicom_load['font'] = helv36
viz_button = tk.Button(button_frame2, text="3D viz", command=perform_3d_visualization)
viz_button['font'] = helv36
eye_segm_button = tk.Button(button_frame3, text="Eye segmentation")
eye_segm_button['font'] = helv36

# Pack buttons into their respective frames
load_button.pack(side=tk.LEFT, padx=10)
dicom_load.pack(side=tk.LEFT, padx=10)
viz_button.pack(padx=10)
eye_segm_button.pack(padx=10)

# Pack frames
button_frame1.pack(pady=10)
button_frame2.pack(pady=10)
button_frame3.pack(pady=10)

visualization_window.mainloop()
