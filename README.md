# apobec
Code used to produce analyses/figures from Buisson et. al. 2019

Main run file is run.github2.m

In order to get the code to run, you must change the path on line 7 to be the full path where the repo was downloaded, e.g.:
srcpath = '/full/path/to/rep/location/';

If you clear the environment variables or exit and resume in the middle of the pipeline, you must reinitialize the srcpath by re-running
lines 7-10 before continuing.

The genome-wide survery of hairpins in lines 255-311 must be edited to ensure the matlab compiler/UGER submission options are compatible
with your system. The resultant files are included in this repo, so this step can be skipped to avoid compatability issues.

