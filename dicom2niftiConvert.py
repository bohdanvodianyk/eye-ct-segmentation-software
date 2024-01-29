import dicom2nifti
import os

patients = os.listdir(f"data/")
for p in patients:
    files = os.listdir(f"data/{p}/")

    for f in files:
        path_dicom = f"data/{p}/{f}"
        print(path_dicom)
        path_nifti = f"output/{p}"
        isExist = os.path.exists(path_nifti)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path_nifti)
            print("The new directory is created!")
        test = dicom2nifti.convert_directory(path_dicom, path_nifti, compression=True)
        print("Done")
