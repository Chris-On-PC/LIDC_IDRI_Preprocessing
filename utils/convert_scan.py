
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


from matplotlib import pyplot as plt


def convert_scan(nods, annotations, scan_id, output_folder, visualize = False, ):

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
