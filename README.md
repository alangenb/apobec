# apobec
Code used to produce analyses/figures from Buisson et. al. 2019

Main run file is run.m

NOTE: This code includes optional sections for compiling + running in parallel to save time.

Example data files are provided in src/data/. This includes a subset of the patient cohort used in the paper that 
has been publicly released previously. Other similarly formatted mutation sets can be used in the same manner, and the input steps in
lines 12-19 should be updated to reflect such changes.  Note, the current code hard-codes the names of the mutational signatures discovered.

Please change path on line 7 to be the full path where the GitHub repo was downloaded to, e.g.:
srcpath = '/full/path/to/rep/location/';

If you clear the environment variables or exit and resume in the middle of the pipeline, you must reinitialize the srcpath by re-running
lines 7-10 before continuing.


