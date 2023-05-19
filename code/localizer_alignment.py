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


def calc_center_of_mass(img_data, affine):
    '''Calc center of mass in scanner coordinates'''
    
    com_ijk = ndimage.center_of_mass(img_data)
    com_ijk = np.array([com_ijk[0], com_ijk[1], com_ijk[2], 1])
    com_xyz = np.matmul(affine, com_ijk)
    return com_xyz

def calc_affine(original_affine, trans_x, trans_y, trans_z, rot_x, rot_y, rot_z, center_of_mass = None):
    '''Function to come up with rotation/translation matrix
    
    Parameters
    ----------
    original_affine : 4x4 numpy array
        Affine matrix of the original image
    trans_x : float
        Translation in x direction
    trans_y : float
        Translation in y direction
    trans_z : float
        Translation in z direction
    rot_x : float
        Rotation in x direction (radians)
    rot_y : float
        Rotation in y direction (radians)
    rot_z : float
        Rotation in z direction (radians)
    center_of_mass : None or 4x1 numpy array
        If None, rotations will be calculated about the axis 0,0,0,
        otherwise of center of mass is provided, the rotations will
        be about this coordinate. The fourth value here should just
        be 1.

    Returns
    -------

    new_affine : 4x4 numpy array
        The new affine after applying the translations and rotations
    
    
    '''
    
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
        flattened_inds_vals = grab_orig_inds_xyz_mat(temp_loc_data, temp_loc_img.affine)
        reshaped_inds = flattened_inds_vals[0].reshape((4, temp_loc_data.shape[0], temp_loc_data.shape[1], temp_loc_data.shape[2]))
        for temp_slice_num in range(num_slices):
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
            
            #Calculate the correlation at the specific slice between the 
            #localizer and full 3d image
            slice_corr = np.round(np.corrcoef(out_vals[~np.isnan(out_vals)].flatten(), temp_slice_data[~np.isnan(out_vals)].flatten())[0,1], 2)
            slice_specific_corrs.append(slice_corr)

            #Plot 
            plt.figure(dpi=100)
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
            plt.imshow(temp_slice_data - out_vals)
            plt.title('Difference (slice corr: {})'.format(slice_corr))
            plt.xticks(np.arange(0,out_vals.shape[0],25), labels='')
            plt.yticks(np.arange(0,out_vals.shape[1],25), labels='')
            plt.gca().grid(color='red', linestyle='-.', linewidth=1)
            plt.tight_layout()
            plt.savefig(os.path.join(output_figures_folder, 'localizer_{}_slice_{}.png'.format(i, temp_slice_num)), bbox_inches='tight')
            #plt.close()
            
    return slice_specific_corrs

def calc_loss(af_vals, localizer_imgs, localizer_vals, reference_data, reference_affine, center_of_mass):
    
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
    num_bins = int((np.max(localizer_vals) - np.min(localizer_vals))/bin_widths)
    bins = np.histogram(localizer_vals, num_bins)
    binned_data = np.zeros(localizer_vals.shape)
    for i in range(bins[0].shape[0] - 1): 
        binned_data[(localizer_vals > bins[0][i+1])*(localizer_vals <= bins[0][i])] = i
    binned_data[localizer_vals >= bins[0][0]] = bins[0].shape[0]
        
    return binned_data

def calc_corr_ratio_loss(af_vals, localizer_imgs, localizer_vals, reference_data, reference_affine, mask_data, center_of_mass, make_plot = False):
    
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
        pass
        
    
    return loss

def make_corr_ratio_loss_plot():
    
    return




def localizer_alignment_anat_update_osprey(anat_files_dict, derivs_folder_path, localizers):



    return