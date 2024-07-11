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

How OSPREY_BIDS is run, and the naming convention for magnetic resonance spectroscopy
files will be determined through the json_settings file that is passed to OSPREY_BIDS
as a positional argument. An example of this type of JSON can be seen in the
codeblock later in this document.

As explained in the XXX section, the json file is a heirarchical dictionary. At
the first level of the heirarchy, the keys describe one grouping of processing
settings. For HBCD these two keys are "HERCULES" and "unedited". These are descriptive
names provided by the user that should be related to the acquisition/processing
details. Further, these names will be propogated to determine the names of the
output folders created by OSPREY_BIDS.

Within each json file there will be a "prerequisites" section. This section can
have two fields including "files" and "files_ref". Each of these fields then has
an associated pattern that will be used with the python glob package to determine
whether a given file in the mrs section of a BIDS dataset matches what is expected
for processing.

.. code-block:: json

   {
    "HERCULES": {
        "seqType": "HERCULES",
        "editTarget": ["GABA","GSH"],
        "MM3coModel": "freeGauss",
        "ECCmetab": "1",
        "prerequisites": {
            "files" : "*_acq-hercules_*svs.nii*",
            "files_ref" : "*_acq-hercules_*ref.nii*"
        }
    },
    "unedited": {
        "seqType": "unedited",
        "ECCmetab": "1",
        "prerequisites": {
            "files" : "*_acq-shortTE_*svs.nii*",
            "files_ref" : "*_acq-shortTE_*ref.nii*"
            }
      }
   }

With the configuration file listed above, if we had the following two
files in the BIDS directory, they would be selected for "HERCULES" processing: ::

   bids_dir/sub-<label>[/ses-<label>]/mrs/sub-<label>[_ses-<label>]_acq-hercules[_run-<label>]_svs.nii.gz
   bids_dir/sub-<label>[/ses-<label>]/mrs/sub-<label>[_ses-<label>]_acq-hercules[_run-<label>]_ref.nii.gz

Further, if we had the following two files in the BIDS directory, they would
be sleected for "unedited" processing: ::

   bids_dir/sub-<label>[/ses-<label>]/mrs/sub-<label>[_ses-<label>]_acq-shortTE[_run-<label>]_svs.nii.gz
   bids_dir/sub-<label>[/ses-<label>]/mrs/sub-<label>[_ses-<label>]_acq-shortTE[_run-<label>]_ref.nii.gz

In both cases, if only one of the two files was present then processing would 
not occur for the configuration missing a file. However, as long as at least one
of the two configurations has all prerequisite files, OSPREY_BIDS processing will be
attempted.

Optional Inputs
===============

Segmentation Directory
----------------------

Beyond the required inputs above, there are two other inputs that can be provided.
These inputs are always used in the Healthy Brain and Child Developement (HBCD) study,
but are not strictly required for processing.

The first optional input is a segmentation directory which is passed to the application
through the "--segmentation_dir" flag. The segmentation directory acts as a BIDS derivatives
directory that has the following files for each subject and session to be processed: ::

   /segmentation_dir/sub-<label>[/ses-<label>]/anat/sub-<label>[_ses-<label>]_space-<modality>_desc-aseg_dseg.nii.gz
   /segmentation_dir/sub-<label>[/ses-<label>]/anat/sub-<label>[_ses-<label>]_space-<modality>_desc-aseg_dseg.json
   /segmentation_dir/sub-<label>[/ses-<label>]/anat/sub-<label>[_ses-<label>]_space-<modality>_desc-brain_mask.nii.gz
   /segmentation_dir/sub-<label>[/ses-<label>]/anat/sub-<label>[_ses-<label>]_space-<modality>_desc-brain_mask.json

Where <modality> is either T1w or T2w, and the nifti images have the same orientation and dimensions
as a T1w or T2w image in the raw BIDS dataset. If the jsons listed above have a SpatialReference field,
that field will be checked to ensure that the T1w (or T2w) image being used for processing has the same
file name as what is listed in the SpatialReference field.

If these segmentation files are present, they will be used to correct for partial tissue effects
within the MRS voxel. To avoid ambiguity in processing outputs, any processing that is attempted
with the "--segmentation_dir" flag will only complete if a segmentation can be identified to use
during processing.

The "aseg" nifti segmentation files listed above should have the same numeric coding as FreeSurfer
"aseg" style images. OSPREY will then recode regions to be grey matter, white matter, or other. Voxels
labeled as other will be assumed to be CSF during the partial volume correction procedure.

If the "--segmentation_dir" flag is not provided, OSPREY will instead utilize SPM segmentation routines
to come up with estimates for grey matter, white matter, and cerebrospinal fluid.

Localizer
---------

By default, OSPREY will make the assumption that the high resolution anatomical (and by extension any
associated segmentations) are registered to the MRS voxels. If you expect this to not be the case for
your data, you should be collecting a localizer scan prior to the MRS acquisition that can be used to
register the MRS voxel to the high resolution anatomical.

If you have localizer images that could be used for this purpose, you should use the
--localizer_registration flag. This will tell OSPREY that you have a localizer image
that should be used for registration purposes. By default, OSPREY will look for localizer
images at the following path: ::

   bids_dir/sub-<label>[/ses-<label>]/anat/[search_term]

In the above example, search term is set by the --localizer_search_term flag and is
*localizer*.nii. In HBCD, the --localizer_search_term value is *mrsLocAx*.nii*. If a given
session directory has more than one localizer image then....