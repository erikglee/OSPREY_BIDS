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

How OSPREY is run, and the naming convention for magnetic resonance spectroscopy
files will be determined through the json_settings file that is passed to OSPREY
as a positional argument. An example of this type of JSON can be seen in the
codeblock later in this document.

As explained in the XXX section, the json file is a heirarchical dictionary. At
the first level of the heirarchy, the keys describe one grouping of processing
settings. For HBCD these two keys are "HERCULES" and "unedited". These are descriptive
names provided by the user that should be related to the acquisition/processing
details. Further, these names will be propogated to determine the names of the
output folders created by OSPREY.

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
of the two configurations has all prerequisite files, OSPREY processing will be
attempted.


