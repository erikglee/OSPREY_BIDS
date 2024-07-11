.. OSPREY_BIDS documentation master file, created by
   sphinx-quickstart on Wed Jun  5 10:48:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Processing Pipeline Details
===========================

As discussed in the expected inputs section, there a number of files that must
be present, following specific formatting requirements, for processing to occur.

At a high level there must be (1) a high resolution anatomical image, (2) MRS images,
and (3) a JSON file that describes how OSPREY should interact with the MRS images 
for processing purposes. 

Beyond that, it is optional that the user also has (4) a FreeSurfer-like segmentation
that is registered to (1), and (5) a localizer image that is registered to (2).

In this section, we describe how OSPREY utilizes (1-5) to do processing. This section
is not meant to describe how input files should be formatted, but rather how they are
used by the OSPREY_BIDS pipeline.


