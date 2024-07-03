.. OSPREY_BIDS documentation master file, created by
   sphinx-quickstart on Wed Jun  5 10:48:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Expected Inputs
===============


This tool expects at minimum one anatomical reference image
and one (or more) images that can be used for MRS analyses.
The anatomical reference image is expected to be under an anat
folder as follows: ::

   bids_dir/sub-<label>[/ses-<label>]/anat/*T1w.nii.gz
   bids_dir/sub-<label>[/ses-<label>]/anat/*T2w.nii.gz

In this example, both a T1w and T2w image are present. Only one of the
two images will be selected for processing. Which image is selected
will depend on whether --preferred_anat_modality has a value of T1w or T2w.
For whichever preferred anat modality is present, there must be exactly one
anatomcial image for the application to choose from. If you also have a
localizer in the anat directory with a T1w/T2w ending, you can utilize
the --terms_not_allowed_in_anat flag to tell the program about what character
sequences can be used to ensure localizer (or other) scans aren't mistakenly
identified as high resolution anatomical images.


The magnetic resonance spectroscopy scans are expected to be in nifti format under
an "mrs" folder as seen in the following example: ::

   bids_dir/suewrtwerwerewr
   bids_dir/swedwerewrewrewr
