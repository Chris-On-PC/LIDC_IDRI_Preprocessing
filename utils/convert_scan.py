
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
from scipy.interpolate import RegularGridInterpolator

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
        vol = None
        for idx, item in enumerate(nod):
            print(item.id)
            
            annoation_volume = annotations.filter(pl.Annotation.id == item.id).first()

            if get_points:
                # vol, mask, pts = annoation_volume.uniform_cubic_resample(side_length=180, return_irp_pts= get_points )
                vol, mask, pts = annoation_volume.custom_cubic_resample(side_length=180, return_irp_pts= get_points )
                get_points = False
                print(vol.shape)
                vol = vol[:180, :180, 58:122]
                print(vol.shape)
            else:
                mask= annoation_volume.custom_cubic_resample(side_length = 180, irp_pts = pts, return_irp_pts= get_points , resample_vol = False)
                mask = mask[:180, :180, 58:122]
            
            # vol = vol[26:154, 26:154, 58:122]
            # mask = mask[26:154, 26:154, 58:122]
            
            
            
            if visualize:

                fig,ax = plt.subplots(1,2,figsize=(5,3))
                print(f'Display slice {int(vol.shape[2]/2)}')
                ax[0].imshow(vol[:,:, int(vol.shape[2]/2)], cmap=plt.cm.gray)
                ax[0].axis('off')

                ax[1].imshow(mask[:,:, int(mask.shape[2]/2)], cmap=plt.cm.gray)
                ax[1].axis('off')

                plt.tight_layout()
                #plt.savefig("../images/mask_bbox.png", bbox_inches="tight")
                plt.show()

            
            

            if "raw" not in h5f.keys():
                ct_vol_to_save = np.swapaxes(vol, 0, 2)
                h5f.create_dataset('raw', data=ct_vol_to_save,  compression="gzip")

            mask = np.swapaxes(mask, 0, 2)
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
        
            

class Custom_annotation(pl.Annotation):

    def __init__(self):
        pass
        # self.name = nae
    def custom_cubic_resample(self, side_length=None, resample_vol=True,
                               irp_pts=None, return_irp_pts=False,
                               resample_img=True, verbose=True, resolution = (0.5, 0.5, 1)):
        """
        Get the CT value volume and respective boolean mask volume. The 
        volumes are interpolated and resampled to have uniform spacing of 1mm
        along each dimension. The resulting volumes are cubic of the 
        specified `side_length`. Thus, the returned volumes have dimensions,
        `(side_length+1,)*3` (since `side_length` is the spacing).
        TODO
        ----
        It would be nice if this function performed fully general 
        interpolation, i.e., not necessarily uniform spacing and allowing 
        different resample-resolutions along different coordinate axes.
        Parameters
        ----------
        side_length: integer, default=None
            The physical length of each side of the new cubic 
            volume in millimeters. The default, `None`, takes the
            max of the nodule's bounding box dimensions.
            If this parameter is not `None`, then it should be 
            greater than any bounding box dimension. If the specified 
            `side_length` requires a padding which results in an 
            out-of-bounds image index, then the image is padded with 
            the minimum CT image value.
        resample_vol: boolean, default=True
            If False, only the segmentation volume is resampled.
        irp_pts: 3-tuple from meshgrid
            If provided, the volume(s) will be resampled over these interpolation
            points, rather than the automatically calculated points. This allows
            for sampling segmentation volumes over a common coordinate-system.
        return_irp_pts: boolean, default=False
            If True, the interpolation points (ix,iy,iz) at which the volume(s)
            were resampled are returned. These can potentially be provided as
            an argument to `irp_pts` for separate selfotations that refer to the
            same nodule, allowing the segmentation volumes to be resampled in a
            common coordinate-system.
        verbose: boolean, default=True
            Turn the loading statement on / off.
        Return
        ------
        [ct_volume,] mask [, irp_pts]: ndarray, ndarray, list of ndarrays
            `ct_volume` and `mask` are the resampled CT and boolean 
            volumes, respectively. `ct_volume` and `irp_points` are optionally
            returned, depending on which flags are set (see above).
        Example
        -------
        An example::
            import numpy as np
            import matplotlib.pyplot as plt
            import pylidc as pl
            ann = pl.query(pl.Annotation).first()
            # resampled volumes will have uniform side length of 70mm and
            # uniform voxel spacing of 1mm.
            n = 70
            vol,mask = ann.uniform_cubic_resample(n)
            # Setup the plot.
            img = plt.imshow(np.zeros((n+1, n+1)), 
                             vmin=vol.min(), vmax=vol.max(),
                             cmap=plt.cm.gray)
            # View all the resampled image volume slices.
            for i in range(n+1):
                img.set_data(vol[:,:,i] * (mask[:,:,i]*0.6+0.2))
                plt.title("%02d / %02d" % (i+1, n))
                plt.pause(0.1)
        """
        bbox  = self.bbox_matrix()
        bboxd = self.bbox_dims()
        rij   = self.scan.pixel_spacing
        rk    = self.scan.slice_spacing

        imin,imax = bbox[0]
        jmin,jmax = bbox[1]
        kmin,kmax = bbox[2]

        xmin,xmax = imin*rij, imax*rij
        ymin,ymax = jmin*rij, jmax*rij

        zmin = self.scan.slice_zvals[kmin]
        zmax = self.scan.slice_zvals[kmax]
        # print(xmin, xmax)
        # print(ymin, ymax)
        # print(zmin, zmax)
        # { Begin input checks.
        if side_length is None:
            side_length = np.ceil(bboxd.max())
        else:
            if not isinstance(side_length, int):
                raise TypeError('`side_length` must be an integer.')
            if side_length < bboxd.max():
                raise ValueError('`side_length` must be greater\
                                   than any bounding box dimension.')
        side_length = float(side_length)
        # } End input checks.

        # Load the images. Get the z positions.
        images = self.scan.load_all_dicom_images(verbose=verbose)
        img_zs = [float(img.ImagePositionPatient[-1]) for img in images]
        img_zs = np.unique(img_zs)

        # Get the z values of the contours.
        contour_zs = np.unique([c.image_z_position for c in self.contours])

        # Get the indices where the nodule stops and starts
        # with respect to the scan z values.
        #kmin = np.where(zmin == img_zs)[0][0]
        #kmax = np.where(zmax == img_zs)[0][0]

        # Initialize the boolean mask.
        mask = self.boolean_mask()

        ########################################################
        # { Begin interpolation grid creation.
        #   (The points at which the volumes will be resampled.)

        # Compute new interpolation grid points in x.
        center_x = (xmax - xmin)/2 + xmin
        d = (side_length/2)*resolution[0]
        
        xhat, step = np.linspace(center_x-d, center_x+d,
                                 int(side_length)+1, retstep=True)
        print(f"x_step_size: {step}")                         
        # assert abs(step-1) < 1e-5, "New x spacing != 1."

        # Do the same for y.
        center_y = (ymax - ymin)/2 + ymin
        d = (side_length/2)*resolution[1]

        yhat, step = np.linspace(center_y-d, center_y+d,
                                 int(side_length)+1, retstep=True)
        print(f"y_step_size: {step}")
        # assert abs(step-1) < 1e-5, "New y spacing != 1."

        # Do the same for z.
        center_z = (zmax - zmin)/2 + zmin
        d = (side_length/2)*resolution[2]
        zhat, step = np.linspace(center_z-d, center_z+d,
                                 int(side_length)+1, retstep=True)
        print(f"z_step_size: {step}")
        # assert abs(step-1) < 1e-5, "New z pixel spacing != 1."

        # } End interpolation grid creation.
        ########################################################

        ########################################################
        # { Begin grid creation.
        #   (The points at which the volumes are assumed to be sampled.)

        # a[x|y|z], b[x|y|z] are the start / stop indexes for the 
        # (non resample) sample grids along each respective axis.

        # It helps to draw a diagram. For example,
        #
        # *--*--*-- ...
        # x3 x4 x5
        #  *---*---*--- ...
        #  xhat0
        #
        # In this case, `ax` would be chosen to be 3
        # since this is the index directly to the left of 
        # `xhat[0]`. If `xhat[0]` is below any grid point,
        # then `ax` is the minimum possible index, 0. A similar
        # diagram helps with the `bx` index.

        T = np.arange(0, 512)*rij

        if xhat[0] <= 0:
            ax = 0
        else:
            ax = (T < xhat[0]).sum() - 1
        if xhat[-1] >= T[-1]:
            bx = 512
        else:
            bx = 512 - (T > xhat[-1]).sum() + 1

        if yhat[0] <= 0:
            ay = 0
        else:
            ay = (T < yhat[0]).sum() - 1
        if yhat[-1] >= T[-1]:
            by = 512
        else:
            by = 512 - (T > yhat[-1]).sum() + 1

        if zhat[0] <= img_zs[0]:
            az = 0
        else:
            az = (img_zs < zhat[0]).sum() - 1
        if zhat[-1] >= img_zs[-1]:
            bz = len(img_zs)
        else:
            bz = len(img_zs) - (img_zs > zhat[-1]).sum() + 1
        
        # These are the actual grids.
        x = T[ax:bx]
        y = T[ay:by]
        z = img_zs[az:bz]
        # print(x,y,z)
        # } End grid creation.
        ########################################################


        # Create the non-interpolated CT volume.
        if resample_vol:
            ctvol = np.zeros(x.shape+y.shape+z.shape, dtype=np.float64)
            for k in range(z.shape[0]):
                ctvol[:,:,k] = images[k+az].pixel_array[ax:bx, ay:by]

        # We currently only have the boolean mask volume on the domain
        # of the bounding box. Thus, we must "place it" in the appropriately
        # sized volume (i.e., `ctvol.shape`). This is accomplished by
        # padding `mask`.
        padvals = [(imin-ax, bx-1-imax), # The `b` terms have a `+1` offset
                   (jmin-ay, by-1-jmax), # from being an index that is
                   (kmin-az, bz-1-kmax)] # corrected with the `-1` here.
        mask = np.pad(mask, pad_width=padvals,
                      mode='constant', constant_values=False)

        # Obtain minimum image value to use as const for interpolation.
        if resample_vol:
            fillval = min([img.pixel_array.min() for img in images])

        if irp_pts is None:
            ix,iy,iz = np.meshgrid(xhat, yhat, zhat, indexing='ij')
        else:
            ix,iy,iz = irp_pts
        IXYZ = np.c_[ix.flatten(), iy.flatten(), iz.flatten()]

        # Interpolate the nodule CT volume.
        if resample_vol:
            rgi = RegularGridInterpolator(points=(x, y, z), values=ctvol,
                                          bounds_error=False, fill_value=fillval)
            ictvol = rgi(IXYZ).reshape(ix.shape)
        # Interpolate the mask volume.
        rgi = RegularGridInterpolator(points=(x, y, z), values=mask,
                                      bounds_error=False, fill_value=False)
        imask = rgi(IXYZ).reshape(ix.shape) > 0

        if resample_vol:
            if return_irp_pts:
                return ictvol, imask, (ix,iy,iz)
            else:
                return ictvol, imask
        else:
            if return_irp_pts:
                return imask, (ix,iy,iz)
            else:
                return imask

