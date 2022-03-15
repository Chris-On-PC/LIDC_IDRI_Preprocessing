
import pylidc as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import os

import h5py
import SimpleITK as sitk  
import nibabel as nib
from PIL import Image
import csv 

from matplotlib import pyplot as plt

def convert_scan_full(scans, output_folder ):

    for scan in scans:

        nods = scan.cluster_annotations()
        print(scan.id)
        ann = pl.query(pl.Annotation).filter(pl.Annotation.scan_id == scan.id)
        # cs.convert_scan(nods, ann, scan.id, output_folder, visualize=False)
        nifti_image = os.path.join(output_folder, "images", f"lung_{scan.id}_0000.nii")
        nifti_labels = os.path.join(output_folder, "labels", f"lung_{scan.id}.nii")
        

        image = scan.to_volume()
        label = np.zeros_like(image)
        transformed_image = np.swapaxes(np.swapaxes(image, 0, 2), 1,2)
        result_image = sitk.GetImageFromArray(transformed_image)
        result_image.SetSpacing(( scan.pixel_spacing, scan.pixel_spacing, scan.slice_spacing))
        result_image = sitk.DICOMOrient(result_image, 'RPI')
        sitk.WriteImage(result_image, nifti_image)

        for id, nod in enumerate(nods):
            for idx, item in enumerate(nod):
                if idx == 0:
                    annoation_volume = ann.filter(pl.Annotation.id == item.id).first()
                    mask = annoation_volume.boolean_mask()
                    label[annoation_volume.bbox()]  = mask
                else:
                    continue
        transformed_label = np.swapaxes(np.swapaxes(label, 0, 2), 1,2).astype(np.int8)            
        result_label = sitk.GetImageFromArray(transformed_label)
        result_label.SetSpacing(( scan.slice_spacing ,scan.pixel_spacing, scan.pixel_spacing))
        result_label = sitk.DICOMOrient(result_label, 'RPI')
        sitk.WriteImage(result_label, nifti_labels, useCompression=True, compressionLevel=9)
        
        print(nifti_labels)

        # writer = sitk.ImageFileWriter()
        # writer.SetFileName(nifti_image)
        # writer.Execute(result_image)

        # writer = sitk.ImageFileWriter()
        # writer.SetFileName(nifti_labels)
        # writer.Execute(result_label)
        # result_label.SetSpacing(scan.slice_spacing)

def convert_scan_full_meta_data(scans, output_folder ):

    for scan in scans:

        nods = scan.cluster_annotations()
        print(scan.id)
        ann = pl.query(pl.Annotation).filter(pl.Annotation.scan_id == scan.id)
        # cs.convert_scan(nods, ann, scan.id, output_folder, visualize=False)
        nifti_image = os.path.join(output_folder, "images", f"lung_{scan.id}_0000.nii")
        nifti_labels = os.path.join(output_folder, "labels", f"lung_{scan.id}.nii")

        meta_path = os.path.join(output_folder, "meta", f"lung_{scan.id}.csv")
        with open(meta_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['id', 'bbox-x', 'bbox-y', 'bbox-z', 'subtlety',
                'internalStructure',
                'calcification',
                'sphericity',
                'margin',
                'lobulation',
                'spiculation',
                'texture',
                'malignancy'])

            image = scan.to_volume()
            label = np.zeros_like(image)
            transformed_image = np.swapaxes(np.swapaxes(image, 0, 2), 1,2)
            result_image = sitk.GetImageFromArray(transformed_image)
            result_image.SetSpacing(( scan.pixel_spacing, scan.pixel_spacing, scan.slice_spacing))
            result_image.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
            # result_image = sitk.DICOMOrient(result_image, 'RPI')
            sitk.WriteImage(result_image, nifti_image)

            for id, nod in enumerate(nods):
                for idx, item in enumerate(nod):
                    if idx == 0:
                        annoation_volume = ann.filter(pl.Annotation.id == item.id).first()
                        mask = annoation_volume.boolean_mask()
                        label[annoation_volume.bbox()]  = mask
                        
                        writer.writerow([id, f"{annoation_volume.bbox()[0].start}-{annoation_volume.bbox()[0].stop}", 
                        f"{label.shape[1]-annoation_volume.bbox()[1].start}-{label.shape[1] - annoation_volume.bbox()[1].stop}", 
                        f"{label.shape[2]- annoation_volume.bbox()[2].start}-{label.shape[2]-annoation_volume.bbox()[2].stop}",
                        annoation_volume.subtlety,
                        annoation_volume.internalStructure,
                        annoation_volume.calcification,
                        annoation_volume.sphericity,
                        annoation_volume.margin,
                        annoation_volume.lobulation,
                        annoation_volume.spiculation,
                        annoation_volume.texture,
                        annoation_volume.malignancy])
                    else:
                        continue

            transformed_label = np.swapaxes(np.swapaxes(label, 0, 2), 1,2).astype(np.int8)            
            result_label = sitk.GetImageFromArray(transformed_label)
            result_label.SetSpacing(( scan.pixel_spacing, scan.pixel_spacing, scan.slice_spacing))
            # result_label.SetSpacing(( scan.slice_spacing ,scan.pixel_spacing, scan.pixel_spacing))
            result_label.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
            # result_label = sitk.DICOMOrient(result_label, 'RPI')
            sitk.WriteImage(result_label, nifti_labels, useCompression=True, compressionLevel=9)
            
            print(nifti_labels)

        # writer = sitk.ImageFileWriter()
        # writer.SetFileName(nifti_image)
        # writer.Execute(result_image)

        # writer = sitk.ImageFileWriter()
        # writer.SetFileName(nifti_labels)
        # writer.Execute(result_label)
        # result_label.SetSpacing(scan.slice_spacing)

def convert_scan_crop(nods, annotations, scan_id, output_folder, visualize = False, ):

    for id, nod in enumerate(nods):

        get_points = True
        pts = None
        h5_filename = os.path.join(output_folder,  f"case_{scan_id}_nodule_{id+1}.h5")
        print(h5_filename)
        h5f = h5py.File(h5_filename,'w')

        for idx, item in enumerate(nod):
            print(item.id)
            
            annoation_volume = annotations.filter(pl.Annotation.id == item.id).first()
            if get_points:
                vol, mask, pts = annoation_volume.uniform_cubic_resample(180, return_irp_pts= get_points )
                get_points = False
            else:
                vol, mask= annoation_volume.uniform_cubic_resample(180, irp_pts =pts, return_irp_pts= get_points )


            vol = vol[26:154, 26:154, 58:122]
            mask = mask[26:154, 26:154, 58:122]
            print(vol.shape)
            if visualize:

                fig,ax = plt.subplots(1,2,figsize=(5,3))

                ax[0].imshow(vol[:,:,32], cmap=plt.cm.gray)
                ax[0].axis('off')

                ax[1].imshow(mask[:,:,32], cmap=plt.cm.gray)
                ax[1].axis('off')

                plt.tight_layout()
                #plt.savefig("../images/mask_bbox.png", bbox_inches="tight")
                plt.show()

            
            vol = np.swapaxes(vol, 0, 2)
            mask = np.swapaxes(mask, 0, 2)

            if "raw" not in h5f.keys():
                h5f.create_dataset('raw', data=vol,  compression="gzip")
            h5f.create_dataset(f'label{idx+1}', data=mask,  compression="gzip")

        # num_labels = len(h5f.keys()) - 1
        # if num_labels < 4:
        #     num_to_add = 4 - num_labels
        #     for label_idx in range(num_to_add):
        #         h5f.create_dataset(f'label{label_idx + num_labels + 1}', data=np.zeros((64, 128, 128)),  compression="gzip")
        #     print(f"Added {num_to_add} empty labels")

        h5f.close()

def convert_h5_to_dicom(filename):
    f = h5py.File(filename, 'r')
    print(f.keys())
    for item in f.keys():
        if 'label' in item:
            img = nib.Nifti1Image(f[item][:].astype(int), np.eye(4))
        else:
            img = nib.Nifti1Image(f[item][:].astype(float), np.eye(4))
        nib.save(img, os.path.join(filename[:-3]+ item + '_converted.nii.gz'))  
    
    f.close()

def h5_to_2D_slices(path, raw_folder, label_folder):

    for item in os.listdir(path):
        print(item)
        h5f = h5py.File(os.path.join(path, item),'r')

        for slice in range(h5f['raw'].shape[0]):
            labels_per_slice = []
            for key in h5f.keys():
                if 'label' in key:
                    label_slice = h5f[key][slice]
                    # print(label_slice.any())
                    if label_slice.any():
                        labels_per_slice.append(label_slice)

            if len(labels_per_slice) != 0:   
               
                ct_img =  h5f['raw'][slice]   
                max = np.amax(ct_img)
                min = np.amin(ct_img)
                # std = np.std(ct_img)
                ct_img = (ct_img - min)/(max - min)
                ct_slice = Image.fromarray(ct_img)
                plt.imsave(os.path.join(raw_folder, f"{item[:-3]}_Slice_{slice}.png"), ct_slice,  cmap=mpl.cm.gray)
                # ct_slice.save(os.path.join(raw_folder, f"{item[:-3]}_Slice_{slice}.png"))
                for idx, label in enumerate(labels_per_slice):
                    label_data = Image.fromarray(label.astype(int)*255)
                    plt.imsave(os.path.join(label_folder, f"{item[:-3]}_Slice_{slice}_annotation_{idx}.png"), label_data, cmap=mpl.cm.gray)
                    # label_data.save(os.path.join(label_folder, f"{item[:-3]}_Slice_{slice}_annotation_{idx}.png"))
        
            

if __name__ == "__main__":
    convert_scan()
