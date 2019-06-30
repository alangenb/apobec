[_Science_ (2019) 364:eaaw2872](https://science.sciencemag.org/content/364/6447/eaaw2872)

## Passenger hotspot mutations in cancer driven by APOBEC3A and mesoscale genomic features

Rémi Buisson<sup>1,2</sup>, Adam Langenbucher<sup>1</sup>, Danae Bowen<sup>2</sup>, Eugene E. Kwan<sup>1</sup>, Cyril H. Benes<sup>1</sup>, Lee Zou<sup>1,3\*</sup>, and Michael S. Lawrence<sup>1,3,4\*</sup>

<sup>1</sup> Massachusetts General Hospital Cancer Center, Harvard Medical School, Boston, Massachusetts, USA.  
<sup>2</sup> Department of Biological Chemistry, University of California - Irvine, Irvine, California, USA.  
<sup>3</sup> Department of Pathology, Massachusetts General Hospital, Harvard Medical School, Boston, Massachusetts, USA.  
<sup>4</sup> Broad Institute of Harvard and MIT, Cambridge, Massachusetts, USA.  

<sup>\*</sup>Corresponding authors:  
Lee Zou, <lzou1@mgh.harvard.edu>  
Michael S. Lawrence, <mslawrence@mgh.harvard.edu>  
Massachusetts General Hospital Cancer Center, Building 149-7th Floor, 13th Street, Charlestown, MA 02129.

# APOBEC hairpins statistical analysis

**Adam Langenbucher and Michael S. Lawrence**

Code used to produce analyses/figures from [Buisson et. al. _Science_ 2019.](https://science.sciencemag.org/content/364/6447/eaaw2872)


Language used is MATLAB.  Built and tested on MATLAB v9.1.0.441655 (R2016b) but should be largely stable across versions. 

Main analysis code has been collected into a single script:  [main run file](run.m)

Helper functions are located here: [src/](src/)

NOTES:

(1) Example data files are provided in [src/data/](src/data/). This includes a subset of the patient cohort used in the paper that 
has been publicly released previously. Other similarly formatted mutation sets can be used in the same manner, and the input steps in
lines 150-153 should be updated to reflect such changes.  Note, the names of the mutational signatures discovered by NMF in the full dataset were hard-coded.  This step requires manual inspection and naming of the signatures that exist in the dataset.

(2) Please change the path on line 3 to be the full path where the GitHub repo was downloaded to, e.g.

    srcpath = '/full/path/to/rep/location/';
If you clear the environment variables or exit and resume in the middle of the pipeline, you must reinitialize the srcpath by re-running lines 3-6 before continuing.
 
(3) This code includes optional sections for compiling and running in parallel to save time.



[![DOI](https://zenodo.org/badge/173772783.svg)](https://zenodo.org/badge/latestdoi/173772783)  **Zenodo snapshot frozen March 2019**

