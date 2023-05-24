#!/usr/bin/env python3
import nibabel as nib
from nibabel import processing
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
from scipy import optimize
import os, glob
import json
import time
import shutil


def calc_center_of_mass(img_data, affine):
    
    com_ijk = ndimage.center_of_mass(img_data)
    com_ijk = np.array([com_ijk[0], com_ijk[1], com_ijk[2], 1])
    com_xyz = np.matmul(affine, com_ijk)
    return com_xyz

def calc_affine(original_affine, trans_x, trans_y, trans_z, rot_x, rot_y, rot_z, center_of_mass = None):
    
    #WHERE TO BRING IN CENTER OF MASS CALCULATION??? MAYBE THINGS WILL CONVERGE FAST
    #ENOUGH WHERE THIS ISN"T NECESSARY
    
    
    #Make empty matrices for rotations
    mat_rot_x = np.eye(4)
    mat_rot_y = np.eye(4)
    mat_rot_z = np.eye(4)
    
    #Pre apply COM mass so that 
    if type(center_of_mass) == type(None):
        temp_COM_mat = np.eye(4)
    else:
        temp_COM_mat = np.eye(4)
        temp_COM_mat[0,3] = -1*center_of_mass[0]
        temp_COM_mat[1,3] = -1*center_of_mass[1]
        temp_COM_mat[2,3] = -1*center_of_mass[2]
    
    #Define mat for x rotations
    mat_rot_x[1,1] = np.cos(rot_x)
    mat_rot_x[2,2] = np.cos(rot_x)
    mat_rot_x[1,2] = -np.sin(rot_x)
    mat_rot_x[2,1] = np.sin(rot_x)
    
    #Define mat for y rotations
    mat_rot_y[0,0] = np.cos(rot_y)
    mat_rot_y[2,2] = np.cos(rot_y)
    mat_rot_y[2,0] = -np.sin(rot_y)
    mat_rot_y[0,2] = np.sin(rot_y)
    
    #Define mat for z rotations
    mat_rot_y[0,0] = np.cos(rot_z)
    mat_rot_y[1,1] = np.cos(rot_z)
    mat_rot_y[0,1] = -np.sin(rot_z)
    mat_rot_y[1,0] = np.sin(rot_z)
    
    #Apply x, then y, then z rotation then add translation
    new_affine = np.matmul(mat_rot_x, temp_COM_mat)
    new_affine = np.matmul(mat_rot_y, new_affine)
    new_affine = np.matmul(mat_rot_z, new_affine)
    new_affine = np.matmul(np.linalg.inv(temp_COM_mat), new_affine)
    #new_affine = np.matmul(mat_rot_y, mat_rot_x)
    #new_affine = np.matmul(mat_rot_z, new_affine)
    new_affine[0,3] = trans_x
    new_affine[1,3] = trans_y
    new_affine[2,3] = trans_z
    #print(new_affine)
    
    return new_affine

def grab_orig_inds_xyz_mat(image_data, affine):
    
    #This should be finished
    
    inds = np.indices(image_data.shape)
    inds_len = inds.shape[1]*inds.shape[2]*inds.shape[3]
    inds_reshaped = np.reshape(inds, (3, inds_len))
    ones = np.ones((1,inds_len))
    full_inds = np.vstack([inds_reshaped, ones])
    orig_xyz = np.matmul(affine, full_inds)
    orig_vals = image_data.flatten().copy()
    
    return orig_xyz, orig_vals

def get_new_xyzs(transformation, original_xyzs):
    
    return np.matmul(transformation, original_xyzs)
    
    
def grab_image_vals(img_data, img_affine, inds, interp_method = 'linear'):
    #Img_data is a matrix you want to sample from. Inds
    #are the xyz inds to grab. out_vals are the interpolated
    #img values at the inds.
    
    
    i = np.arange(0,img_data.shape[0])
    j = np.arange(0,img_data.shape[1])
    k = np.arange(0,img_data.shape[2])
    interp = RegularGridInterpolator((i, j, k), img_data, method = interp_method, bounds_error = False)
    inds_xyz_to_ijk = np.matmul(np.linalg.inv(img_affine), inds)
    
    
    out_vals = interp(inds_xyz_to_ijk[0:3,:].transpose())
    return out_vals

   
def make_alignment_images(full_registered_nifti_path, localizers_arr_path, output_figures_folder, close_figures = True):
    '''Make plots to show overlap between some full volumetric nifti and localizers
    
    This function makes overlays to show the alignment between the image
    represented by full_registered_nifti_path and the images within
    localizers_arr_path. The overlays are saved to output_figures_folder
    and 
    
    Parameters
    ----------
    
    full_registered_nifti_path : str
        Path to full 3d nifti image, probably an image that
        is registered to the localizers
    localizers_arr_path : list
        List of paths to localizer images that will be used
        to generate overlays
    output_figures_folder : str
        Path to the folder to store overlays. This will be
        created if it doesn't already exist
        
    Returns
    -------
    
    slice_specific_corrs : list
        List of slice specific corrs between the full volumetric
        image and a slice from any of the localizer images
    
    '''
    
    full_img = nib.load(full_registered_nifti_path)
    full_data = full_img.get_fdata()
    full_affine = full_img.affine
    
    i = np.arange(0,full_data.shape[0])
    j = np.arange(0,full_data.shape[1])
    k = np.arange(0,full_data.shape[2])
    interp = RegularGridInterpolator((i, j, k), full_data, method = 'linear', bounds_error = False)
    
    slice_specific_corrs = []
    
    if os.path.exists(output_figures_folder) == False:
        os.makedirs(output_figures_folder)

    for i, temp_localizer in enumerate(localizers_arr_path):

        temp_loc_img = nib.load(temp_localizer)
        temp_loc_data = temp_loc_img.get_fdata()
        smallest_dim = np.argmin(temp_loc_data.shape)
        num_slices = temp_loc_data.shape[smallest_dim]
        if num_slices > 4:
            slices = np.round(np.linspace(0, num_slices, 6)[1:-2]).astype(int)
        else:
            slices = np.linspace(0,num_slices - 1,num_slices).astype(int)
        flattened_inds_vals = grab_orig_inds_xyz_mat(temp_loc_data, temp_loc_img.affine)
        reshaped_inds = flattened_inds_vals[0].reshape((4, temp_loc_data.shape[0], temp_loc_data.shape[1], temp_loc_data.shape[2]))
        for temp_slice_num in slices:
            if smallest_dim == 0:
                temp_slice = reshaped_inds[:,temp_slice_num,...]
                temp_slice_data = temp_loc_data[temp_slice_num,:,:]
            elif smallest_dim == 1:
                temp_slice = reshaped_inds[:,:,temp_slice_num,...]
                temp_slice_data = temp_loc_data[:,temp_slice_num,:]
            elif smallest_dim == 2:
                temp_slice = reshaped_inds[:,:,:,temp_slice_num]
                temp_slice_data = temp_loc_data[:,:,temp_slice_num]
            else:
                raise ValueError('Error: localizer should be 3d image')
            flattened_slice_inds = temp_slice.reshape((temp_slice.shape[0], int(temp_slice.shape[1]*temp_slice.shape[2])))


            #Find values in the full 3d image that correspond to the 
            #current slice in the current localizer image
            inds_xyz_to_ijk = np.matmul(np.linalg.inv(full_affine), flattened_slice_inds)
            out_vals = interp(inds_xyz_to_ijk[0:3,:].transpose())
            out_vals = out_vals.reshape((temp_slice.shape[1], temp_slice.shape[2]))

            #Plot 
            plt.figure(dpi=200)
            plt.subplot(1,3,1)
            plt.imshow(out_vals)
            plt.xticks(np.arange(0,out_vals.shape[0],25), labels='')
            plt.yticks(np.arange(0,out_vals.shape[1],25), labels='')
            plt.gca().grid(color='red', linestyle='-.', linewidth=1)
            plt.title('Full Volumetric Img.')
            plt.subplot(1,3,2)
            plt.imshow(temp_slice_data)
            plt.title('Localizer Img.')
            plt.xticks(np.arange(0,out_vals.shape[0],25), labels='')
            plt.yticks(np.arange(0,out_vals.shape[1],25), labels='')
            plt.gca().grid(color='red', linestyle='-.', linewidth=1)
            plt.subplot(1,3,3)
            plt.imshow((temp_slice_data - np.nanmean(temp_slice_data))/np.nanstd(temp_slice_data) - (out_vals - np.nanmean(out_vals))/np.nanstd(out_vals))
            plt.title('Difference')
            plt.xticks(np.arange(0,out_vals.shape[0],25), labels='')
            plt.yticks(np.arange(0,out_vals.shape[1],25), labels='')
            plt.gca().grid(color='red', linestyle='-.', linewidth=1)
            plt.tight_layout()
            plt.savefig(os.path.join(output_figures_folder, 'localizer_{}_slice_{}.png'.format(i, temp_slice_num)), bbox_inches='tight')
            if close_figures:
                plt.close()                
    return

def calc_loss(af_vals, localizer_imgs, localizer_vals, reference_data, reference_affine, center_of_mass, xyz_s_list):
    
    affine_transforms = []
    new_xyz_s_list = []
    for i, temp_img in enumerate(localizer_imgs):
        affine_transforms.append(calc_affine(localizer_imgs[i].affine, af_vals[0], af_vals[1], af_vals[2], af_vals[3], af_vals[4], af_vals[5], center_of_mass = center_of_mass)) #transform to apply
        new_xyz_s_list.append(get_new_xyzs(affine_transforms[i], xyz_s_list[i]))
        
    new_xyz_s_arr = np.hstack(new_xyz_s_list)
    reference_vals = grab_image_vals(reference_data, reference_affine, new_xyz_s_arr, interp_method = 'linear')
    good_ref = reference_vals[np.isnan(reference_vals) == False]
    good_loc = localizer_vals[np.isnan(reference_vals) == False]
    loss = 1 - np.corrcoef(good_ref, good_loc)[0,1]
    
    return loss

def calc_localizer_val_bins(localizer_vals):
    
    std = np.std(localizer_vals)
    bin_widths = 3.49*std*np.power(localizer_vals.shape[0], -1/3) #This equation is optimal for unimodal case per page 151 of jenkinson paper 
    #bin_widths = 2*std*np.power(localizer_vals.shape[0], -1/3)
    num_bins = int((np.max(localizer_vals) - np.min(localizer_vals))/bin_widths)
    bins = np.histogram(localizer_vals, num_bins)
    binned_data = np.zeros(localizer_vals.shape)
    for i in range(bins[0].shape[0] - 1): 
        binned_data[(localizer_vals > bins[0][i+1])*(localizer_vals <= bins[0][i])] = i
    binned_data[localizer_vals >= bins[0][0]] = bins[0].shape[0]
        
    return binned_data

def calc_corr_ratio_loss(af_vals, localizer_imgs, localizer_vals, reference_data, reference_affine, mask_data, center_of_mass, xyz_s_list, make_plot = False, image_output_path = None):
    
    affine_transforms = []
    new_xyz_s_list = []
    for i, temp_img in enumerate(localizer_imgs):
        affine_transforms.append(calc_affine(localizer_imgs[i].affine, af_vals[0], af_vals[1], af_vals[2], af_vals[3], af_vals[4], af_vals[5], center_of_mass = center_of_mass)) #transform to apply
        new_xyz_s_list.append(get_new_xyzs(affine_transforms[i], xyz_s_list[i]))
        
    new_xyz_s_arr = np.hstack(new_xyz_s_list)
    reference_vals = grab_image_vals(reference_data, reference_affine, new_xyz_s_arr, interp_method = 'linear')
    mask_vals = grab_image_vals(mask_data, reference_affine, new_xyz_s_arr, interp_method = 'nearest')
    good_ref = reference_vals[(np.isnan(reference_vals) == False)*(mask_vals > 0.5)]
    good_loc = localizer_vals[(np.isnan(reference_vals) == False)*(mask_vals > 0.5)]
    num_good = good_loc.shape[0]
    
    unique_loc_vals = np.unique(good_loc)
    corr_ratio = 0
    for i in range(unique_loc_vals.shape[0]):
        n_k = np.sum(unique_loc_vals == unique_loc_vals[i])
        corr_ratio += (n_k/num_good)*np.var(good_ref[good_loc == unique_loc_vals[i]])
    corr_ratio = corr_ratio/np.var(good_ref)
    
    loss = 1 - corr_ratio
    loss = corr_ratio
    
    
    if make_plot:
        make_corr_ratio_loss_plot(unique_loc_vals, good_loc, good_ref, output_image_path = image_output_path, close_image = False)
        
    
    return loss

def make_corr_ratio_loss_plot(unique_loc_vals, good_loc, good_ref, output_image_path = None, close_image = True):
    
    plt.figure(dpi = 100)
    differences_1 = []
    bin_jitters_1 = []
    differences_2 = []
    bin_jitters_2 = []
    for i in range(unique_loc_vals.shape[0]):
        
        temp_vals = good_ref[good_loc == unique_loc_vals[i]]
        temp_differences = np.log10(np.absolute(temp_vals - np.mean(temp_vals)))*np.sign(temp_vals - np.mean(temp_vals))
        temp_bin_jitters = np.random.uniform(low = i, high = i + 1, size = temp_differences.shape)
        
        if np.mod(i,2) == 0:
            differences_1.append(temp_differences)
            bin_jitters_1.append(temp_bin_jitters)
        else:
            differences_2.append(temp_differences)
            bin_jitters_2.append(temp_bin_jitters)
    
    differences_1 = np.hstack(differences_1)
    bin_jitters_1 = np.hstack(bin_jitters_1)
    differences_2 = np.hstack(differences_2)
    bin_jitters_2 = np.hstack(bin_jitters_2)
    plt.scatter(bin_jitters_1, differences_1, s = 0.05)
    plt.scatter(bin_jitters_2, differences_2, s = 0.05)
    plt.xlabel('Bin Number')
    plt.ylabel('sign(Deviation)*log10(absolute(Deviation))')
    plt.axhline(2.5, linestyle = '--', color = 'grey', linewidth = 1)
    plt.axhline(-2.5, linestyle = '--', color = 'grey', linewidth = 1)
    if type(output_image_path) == type(None):
        pass
    else:
        plt.savefig(output_image_path)
    
    if close_image:
        plt.close()
        
    return

def make_readme(affine_readme_path):
    
    with open(affine_readme_path, 'w') as f:
        f.write('This folder contains registration results from a high res anatomical template to a localizer thats presumed to be in MRS voxel space.\n')
        f.write('The details of what images were registered can be found in registration_summary.json and figures showing the quality of the registration can be found in the figures folder.\n')
        f.write('The new copy of the reference image (now aligned to the localizer) is found at reference_img_aligned_to_localizer.nii.gz\n\n')
        f.write('How to use transform_mat.npy file:\n\n\n')
        f.write('import nibabel as nib\nimport numpy as np\npath_to_image_in_reference_space = ""\ntemp_img = nib.load(path_to_image_in_reference_space)\n')
        f.write('transform_mat = np.load("transform_mat.npy")\ntemp_img.affine = np.matmul(transform_mat, temp_img.affine)\n')
        f.write('nib.save(temp_img, "/some/new/path/for/image/now/in/localizer/space/img.nii.gz")')
        
    return


def localizer_alignment_anat_update_osprey(anat_files_dict, registration_output_folder, localizer_paths):
    '''Registers anat reference image to the localizer image(s)
    
    
    
    Parameters
    ----------
    
    anat_files_dict : dict
        Has (at minimum) key 'files_nii' and optionally 'files_seg' that will be
        registered to the localizer image. After registration,
        this path will be reset to be the path to the new image
        following registration.
    registration_output_folder : str
        Path to the folder that will be created to store registration
        results. This will be subject/ses and in certain cases run
        specific.
    localizer_paths : list of strings
        The paths to the localizer image or images to be registered.
        You will only have multiple entries in this list if axial
        images were stored in different images than sagital or coronal.
        
    Returns
    -------
    
    anat_files_dict : dict
        The same dictionary as before, but now the 'files_nii' key
        has been updated to point to the registered image
    
    
    
    '''
    
    output_folder = registration_output_folder
    if os.path.exists(os.path.join(output_folder, 'figures')) == False:
        os.makedirs(os.path.join(output_folder, 'figures'))
    make_readme(os.path.join(output_folder, 'readme.txt'))

    #Load the reference image
    reference_path = anat_files_dict['files_nii'][0]
    reference_img = nib.load(reference_path)
    reference_data = reference_img.get_fdata()

    #Load and dilate brain mask by 10 iterations ... this keeps the scalp in registration but not neck
    #USING MASK SEEMED TO HURT REGISTRATION FOR LOWRES LOCALIZERS SO AM EXCLUDING THIS FOR NOW
    #mask_data = nib.load(brain_mask_path).get_fdata()
    #mask_data = ndimage.binary_dilation(mask_data, iterations = 10)
    #mask_data = mask_data.astype(float) + 1

    reference_data_10mm_smoothing = processing.smooth_image(reference_img, 10).get_fdata()
    reference_com = calc_center_of_mass(reference_data, reference_img.affine)
    mask_data = np.ones(reference_data.shape) #we arent using a mask right now so this is just a dummy mask
    #reference_com = None

    localizer_imgs = []
    xyz_s_list = []
    vals = []
    localizer_sizes = []
    for i, temp_path in enumerate(localizer_paths):
        localizer_imgs.append(nib.load(temp_path))
        temp_xyz, temp_vals = grab_orig_inds_xyz_mat(localizer_imgs[i].get_fdata(), localizer_imgs[i].affine)
        xyz_s_list.append(temp_xyz)
        vals.append(temp_vals)

        localizer_sizes.append(localizer_imgs[i].get_fdata().size)

    xyz_s_arr = np.hstack(xyz_s_list)
    localizer_vals = np.hstack(vals)
    localizer_vals = calc_localizer_val_bins(localizer_vals) #NOW THESE ARE BINS
    reference_vals = grab_image_vals(reference_data, reference_img.affine, xyz_s_arr, interp_method = 'linear')


    good_ref = reference_vals[np.isnan(reference_vals) == False]
    good_loc = localizer_vals[np.isnan(reference_vals) == False]
    print('Original Ref/Localizer Correlation Ratio (0 is best, 1 is worst):')
    original_corr = calc_corr_ratio_loss([0,0,0,0,0,0], localizer_imgs, localizer_vals, reference_data, reference_img.affine, mask_data, reference_com, xyz_s_list, make_plot = True, image_output_path = os.path.join(registration_output_folder, 'figures', 'corr_ratio_pre_registration.png'))
    print(original_corr)

    bounds_10mm = [[-100,100],[-100,100],[-100,100],[-1.5,1.5],[-1.5,1.5],[-1.5,1.5]]
    tic = time.perf_counter()
    options = {'maxfun':5000, 'maxiter':50}
    results_tnc_10mm = optimize.minimize(calc_corr_ratio_loss, [0,0,0,0,0,0], args=(localizer_imgs, localizer_vals, reference_data_10mm_smoothing, reference_img.affine, mask_data, reference_com, xyz_s_list),
                                          method='TNC', jac=None, bounds=bounds_10mm, options=options)
    results_tnc_00mm = optimize.minimize(calc_corr_ratio_loss, results_tnc_10mm.x, args=(localizer_imgs, localizer_vals, reference_data, reference_img.affine, mask_data, reference_com, xyz_s_list),
                                          method='TNC', jac=None, bounds=bounds_10mm, options=options)

    toc = time.perf_counter()
    print(f"Ran optimization in {toc - tic:0.4f} seconds")

    ###Illustrate the performance of the new transformation
    affine_transforms = []
    new_xyz_s_list = []
    for i, temp_img in enumerate(localizer_imgs):
        affine_transforms.append(calc_affine(localizer_imgs[i].affine, results_tnc_00mm.x[0], results_tnc_00mm.x[1], results_tnc_00mm.x[2], results_tnc_00mm.x[3], results_tnc_00mm.x[4], results_tnc_00mm.x[5], reference_com)) #transform to apply
        new_xyz_s_list.append(get_new_xyzs(affine_transforms[i], xyz_s_list[i]))
    new_xyz_s_arr = np.hstack(new_xyz_s_list)
    reference_vals = grab_image_vals(reference_data, reference_img.affine, new_xyz_s_arr, interp_method = 'linear')
    good_ref = reference_vals[np.isnan(reference_vals) == False]
    good_loc = localizer_vals[np.isnan(reference_vals) == False]
    registered_corr = calc_corr_ratio_loss(results_tnc_00mm.x, localizer_imgs, localizer_vals, reference_data, reference_img.affine, mask_data, reference_com, xyz_s_list, make_plot = True, image_output_path = os.path.join(registration_output_folder, 'figures', 'corr_ratio_post_registration.png'))
    if original_corr > registered_corr:
        inv_affine = calc_affine(np.eye(4), results_tnc_00mm.x[0], results_tnc_00mm.x[1], results_tnc_00mm.x[2], results_tnc_00mm.x[3], results_tnc_00mm.x[4], results_tnc_00mm.x[5], reference_com)
        inv_affine = np.linalg.inv(inv_affine)
    else:
        inv_affine = np.eye(4)
        registered_corr = original_corr
        shutil.copyfile(os.path.join(registration_output_folder, 'figures', 'corr_ratio_pre_registration.png'), os.path.join(registration_output_folder, 'figures', 'corr_ratio_post_registration.png'))
        
    print('Registered Ref/Localizer Correlation (0 is best, 1 is worst):')
    print(registered_corr)


    if original_corr > registered_corr:
        inv_affine = calc_affine(np.eye(4), results_tnc_00mm.x[0], results_tnc_00mm.x[1], results_tnc_00mm.x[2], results_tnc_00mm.x[3], results_tnc_00mm.x[4], results_tnc_00mm.x[5], reference_com)
        inv_affine = np.linalg.inv(inv_affine)
    else:
        inv_affine = np.eye(4)
        registered_corr = original_corr
    new_affine = np.matmul(inv_affine, reference_img.affine)
    new_img = nib.nifti1.Nifti1Image(reference_data, new_affine)
    registered_output_image_name = os.path.join(output_folder, 'reference_img_aligned_to_localizer.nii.gz')
    nib.save(new_img, registered_output_image_name)
    np.save(os.path.join(output_folder, 'transform_mat.npy'), inv_affine) #This can be used with other images to update their affines
    if 'files_seg' in anat_files_dict.keys():
        new_seg_img = nib.load(anat_files_dict['files_seg'])
        new_seg_img = nib.Nifti1Image(new_seg_img.get_fdata(), new_affine)
        print('Saving new segmentation image.')
        nib.save(new_seg_img, os.path.join(output_folder, 'reference_seg_aligned_to_localizer.nii.gz'))
        anat_files_dict['files_seg'] = os.path.join(output_folder, 'reference_seg_aligned_to_localizer.nii.gz')


    make_alignment_images(registered_output_image_name, localizer_paths, os.path.join(output_folder, 'figures'))
    registration_dict = {"reference_img" : reference_path,
                         "localizer_imgs": localizer_paths,
                         "reference_localizer_corr_ratio_pre_registration" : np.round(original_corr, 8),
                         "reference_localizer_corr_ratio_post_registration" : np.round(registered_corr, 8)}

    with open(os.path.join(output_folder, 'registration_summary.json'), 'w') as f:
        f.write(json.dumps(registration_dict, indent = 6))
        
    anat_files_dict['files_nii'] = registered_output_image_name
        
    return anat_files_dict