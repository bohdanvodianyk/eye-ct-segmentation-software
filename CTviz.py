import tkinter as tk
from tkinter import filedialog
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def apply_threshold_and_normalize(slice, lower_threshold, upper_threshold):
    # Apply threshold: set values outside the threshold range to 0
    thresholded_slice = np.copy(slice)
    thresholded_slice[slice < lower_threshold] = 0
    thresholded_slice[slice > upper_threshold] = 0

    # Normalize the remaining values
    valid_mask = (thresholded_slice > 0)
    if np.any(valid_mask):
        normalized_slice = np.zeros_like(thresholded_slice, dtype=np.float32)
        normalized_slice[valid_mask] = (thresholded_slice[valid_mask] - np.min(slice)) / (
                    np.max(slice) - np.min(slice))
        normalized_slice *= 255
        return normalized_slice
    else:
        return thresholded_slice


def load_and_visualize_nifti():
    global file_path, image_array, lower_threshold, upper_threshold

    file_path = filedialog.askopenfilename(filetypes=[("NIfTI files", "*.nii;*.nii.gz")])
    if not file_path:
        return

    nifti_image = nib.load(file_path)
    image_array = nifti_image.get_fdata()

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.35)

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
        slice_idx = int(slice_slider.val)
        lower_threshold = lower_slider.val
        upper_threshold = upper_slider.val
        normalized_image = apply_threshold_and_normalize(image_array[:, :, slice_idx], lower_threshold, upper_threshold)
        img.set_data(normalized_image)
        fig.canvas.draw_idle()

    slice_slider.on_changed(update)
    lower_slider.on_changed(update)
    upper_slider.on_changed(update)

    plt.show()

def perform_3d_visualization():
    global image_array, lower_threshold, upper_threshold
    if image_array is not None:
        # Implement 3D visualization using VTK here
        print("3D visualization to be implemented")

def open_3d_visualization_window():
    visualization_window = tk.Toplevel()
    visualization_window.title("3D Visualization")

    load_button = tk.Button(visualization_window, text="Load NIfTI File", command=load_and_visualize_nifti)
    load_button.pack(pady=10)

    viz_button = tk.Button(visualization_window, text="3D viz", command=perform_3d_visualization)
    viz_button.pack(pady=10)

root = tk.Tk()
root.title("NIfTI Viewer")

open_button = tk.Button(root, text="Open NIfTI File", command=load_and_visualize_nifti)
open_button.pack(pady=20)

root.mainloop()
