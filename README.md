# CrispressoTools Ver 1.0.1

Tools to facilitate Crispresso Analysis. Fixes and updates from `as114`

## fixNonOverlappingReads: 
Script to allow non-overlapping paired end reads to
work with Crispresso's paired end mode.

1. Added logic to extend reads to a length that is calculated from length of reference amplicon sequence, so that the maximum overlap between read1 and read2 is 10bp. Crispresso software penalizes reads if overlap is more than certain number of base pairs.

2. Added new script CrispressoPipeline.py to submit jobs to Grid engine. This is the script that will be used to start CRISPR analysis.

3. Commented the code for clarity.

4. The code has been tested with several projects and the results were compared with results from old script without padding.
