# APOBEC hairpins

Passenger hotspot mutations in cancer driven by APOBEC3A and mesoscale genomic features

Rémi Buisson1,2, Adam Langenbucher1, Danae Bowen2, Eugene E. Kwan1, Cyril H. Benes1, 
Lee Zou1,3\*, and Michael S. Lawrence1,3,4\*

1) Massachusetts General Hospital Cancer Center, Harvard Medical School, Boston, Massachusetts, USA.
2) Department of Biological Chemistry, University of California, Irvine, California, USA. 
3) Department of Pathology, Massachusetts General Hospital, Harvard Medical School, Boston, Massachusetts, USA.
4) Broad Institute of Harvard and MIT, Cambridge, Massachusetts, USA.

\*Corresponding authors:

Lee Zou, Massachusetts General Hospital Cancer Center, Building 149-7th Floor, 13th Street, Charlestown, MA 02129. Phone: 617-724-9534; Fax: 617-726-7808; E-mail: lzou1@mgh.harvard.edu

Michael S. Lawrence, Massachusetts General Hospital Cancer Center, Building 149-7th Floor, 13th Street, Charlestown, MA 02129. Phone 617-643-4379; Fax: 617-726-7808; E-mail: lawrence@broadinstitute.org    

=====================================================================


# APOBEC hairpins statistical analysis

Adam Langenbucher and Mike Lawrence

Code used to produce analyses/figures from Buisson et. al. 2019

Language used is MATLAB.  Built and tested on MATLAB v9.1.0.441655 (R2016b) but should be largely stable across versions. 

All analysis code has been collected into a single script:

[main run file](run.m)

NOTES:

(1) Example data files are provided in src/data/. This includes a subset of the patient cohort used in the paper that 
has been publicly released previously. Other similarly formatted mutation sets can be used in the same manner, and the input steps in
lines 150-153 should be updated to reflect such changes.  Note, the names of the mutational signatures that exist in the sample dataset have been hard-coded here.  When analyzing any other dataset, this step will require manual inspection and naming of the signatures that exist in that dataset.

(2) Please change the path on line 3 to be the full path where the GitHub repo was downloaded to, e.g.

    srcpath = '/full/path/to/rep/location/';
If you clear the environment variables or exit and resume in the middle of the pipeline, you must reinitialize the srcpath by re-running lines 3-6 before continuing.
 
(3) This code includes optional sections for compiling + running in parallel to save time.


# Zenodo snapshot frozen March 2019
[![DOI](https://zenodo.org/badge/173772783.svg)](https://zenodo.org/badge/latestdoi/173772783)

