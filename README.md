# apobec
Code used to produce analyses/figures from Buisson et. al. 2019

Main run file is run.github2.m

Example data files are provided in src/data/. This includes a subset of the patient cohort used in the paper that 
is available for public release. Other similarly formatted mutation sets can be used in the same manner, and the input steps in
lines 12-19 should be updated to reflect such changes.

In order to get the code to run, you must change the path on line 7 to be the full path where the repo was downloaded, e.g.:
srcpath = '/full/path/to/rep/location/';

If you clear the environment variables or exit and resume in the middle of the pipeline, you must reinitialize the srcpath by re-running
lines 7-10 before continuing.

The genome-wide survery of hairpins in lines 255-307 must be edited to ensure the matlab compiler/UGER submission options are compatible
with your system. The resultant files are included in this repo, so this step can be skipped to avoid compatability issues.

