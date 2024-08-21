.. OSPREY_BIDS documentation master file, created by
   sphinx-quickstart on Wed Jun  5 10:48:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Default processing options
==========================

The minimum input parameters for the Osprey analysis are shown in the JSON Configuration section.
For "HERCULES" scans, the parameters "seqType" (default value = "HERCULES"),
"editTarget" (default value = ["GABA", "GSH"]), and "ECCmetab" (default value = "1") should
not be modified as they are required for correct processing.
Similarly, "seqType" (default value = "UnEdited") and "ECCmetab" (default value = "1")
should also not be modified for the "unedited" short-TE datasets. All parameter options
that users can set in Osprey job files are explained in detail in the documentation (https://schorschinho.github.io/osprey/)
and example job files (in the ‘exampledata’ folder at https://github.com/schorschinho/osprey)
and in the original Osprey publication (https://doi.org/10.1016/j.jneumeth.2020.108827).

Additional processing options
==========================

There are several key parameters that can be modified to change different aspects (pre-processing and modeling)
of the analysis workflow. The following section contains a brief description of the available options:

Spectral registration
---------------------

This option sets the algorithm Osprey will use to align individual spectral transients for improved linewidth and signal-to-noise ratio. ::

	"SpecReg": "RobSpecReg"
	"SpecReg": "ProbSpecReg"
	"SpecReg": "RestrSpecReg"
	"SpecReg": "none"

"RobSpecReg" (default) is a spectral alignment of the transients with water/lipid removal,
similarity metric, and weighted averaging (https://doi.org/10.1002/nbm.4368). "ProbSpecReg"
is probabilistic spectral alignment to median target and weighted averaging (https://doi.org/10.1002/mrm.27027).
"RestrSpecReg" is frequency-restricted (to fit range) spectral alignment using similarity metric and weighted averaging.
"None" uses no spectral alignment of the transients.

Subspectra alignment
--------------------

This option is for "HERCULES" only. It determines how Osprey aligns the four sub-experiments of the HERCULES acquisition to each other. ::

	"SubSpecAlignment": "L2Norm"
	"SubSpecAlignment": "L1Norm"
	"SubSpecAlignment": "none"


"L2Norm" (default) performs optimization using L2-norm to minimize the difference between
different pairs of sub-spectra over various frequency ranges for optimal alignment. "L1Norm"
uses the L1-norm to minimize the difference between pairs of sub-spectra over the range between
1.95 and 4 ppm (https://doi.org/10.1016/j.jmr.2017.04.004). "none" applies no sub-spectra alignment. ::

	"UnstableWater": "0"
	"UnstableWater": "1"

If set to "0"  (default), the default L2Norm sub-spectra alignment (see above) is performed using
the residual water peak (~4.68 ppm) to align sub-spectra A and B. If set to "1", L2Norm sub-spectra
alignment is instead performed using the total choline peak at 3.22 ppm. This can improve their
alignment in the case of unstable residual water.

Fit options
-----------

There are several options to change the modeling behavior.
It is best to keep those at the default settings since MRS
fit settings are well-known to substantially affect the quantitative
results. However, they are available for advanced MRS researchers to modify. 

Several model setting parameters describe the co-edited macromolecule (MM)
model at 3 ppm in the HERCULES GABA-edited difference spectrum (https://doi.org/10.1002/nbm.4618): ::


	"MM3coModel": "freeGauss"
	"MM3coModel": "fixedGauss"
	"MM3coModel": "3to2MM"
	"MM3coModel": "3to2MMsoft"
	"MM3coModel": "1to1GABA"
	"MM3coModel": "1to1GABAsoft"
	"MM3coModel": "none"

- "freeGauss" (default) defines the 3-ppm MM peak in the GABA-edited difference
  spectrum as a 2-proton Gaussian peak with a free full-width at half-maximum (FWHM)
  parameter to be estimated during modeling.
- "fixedGauss" defines the 3-ppm MM peak in the GABA-edited difference spectrum as a 
  2-proton Gaussian peak with a fixed FWHM of 14 Hz.
- "3to2MM" defines a composite MM basis function combining the 3-ppm and the 0.9-ppm MM
  peaks with an amplitude ratio of 3 to 2. 
- "3to2MMsoft" employs soft constraints during modeling to enforce an amplitude ratio of 3 to 2
  for the 3-ppm and the 0.9-ppm MM peak.
- "1to1GABA" defines a composite basis function combining the 3-ppm MM peak and the GABA
  basis function. "1to1GABAsoft" employs soft constraints during modeling to enforce an
  amplitude ratio of 1 to 1 for the 3-ppm MM peak and the GABA basis function.
  More information can be found here (https://doi.org/10.1002/nbm.4618)

The next option allows to change the FWHM (in Hz) of the Gaussian peak of the 3-ppm MM peak: ::

	"FWHMM3co": "14"

Finally, a set of options governs the spectral fit range
(the frequency over which the model is optimized) as well
as the stiffness of the baseline.

The option to modify the metabolite model frequency range in ppm: ::


	"lolim_range": "0.5",
	"uplim_range": "4.2"

The option to modify the water model frequency range in ppm: ::

	"lolim_rangew": "2.0",
	"uplim_rangew": "7.4"


The option to modify the minimal spacing of neighboring knots of the cubic spline baseline: ::

	"bLineKnotSpace": "0.4"


Output files
============

Osprey generates several derivative files. The most interesting analysis results,
namely the metabolite estimates, can be found in the `QuantifyResults` folder.
This folder contains tab-separated value (.tsv) files with the analysis results using
different quantification methods. For HERCULES, these files are generated for each
modeled sub-spectrum (diff1, diff2, sum). Each .tsv file is accompanied by a matching
.json file which holds more detailed explanations of the exact quantification process.
Please consult the original Osprey manuscript for further details. Additional quality
metrics (linewidth, signal-to-noise ratio, etc.) can be found in the QM_processed_spectra.tsv
file, again accompanied by a .json descriptor file.