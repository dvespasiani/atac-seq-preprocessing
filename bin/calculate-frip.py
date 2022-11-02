## Author: Laura E Cook, University of Melbourne
## Purpose:
#   1. Intersect aligned reads with called peaks
#   2. Count number of lines in intersected file
#   3. Count number of lines in aligned reads (BED) file
#   4. Calculate fraction of peaks that are in the aligned reads
#   5. ENCODE guidelines state that FRiP > 0.3 is optimal, FRiP > 0.2 is acceptable

import sys
import os
sys.path.append('/home/users/allstaff/vespasiani.d/.local/lib/python3.6/site-packages')
from pybedtools import BedTool

# bed file (converted from aligned BAM file to bed)
bed = BedTool(sys.argv[1])

# called narrowPeak file from MACS2
peak = BedTool(sys.argv[2])

# overlap bed and peak files
overlap = bed.intersect(peak, nonamecheck=True, wa=True, u=True)

num_lines_overlap = 0

# determine how many lines in overlap file
for line in overlap:
    num_lines_overlap += 1
print ('Number of reads that intersect with peaks = ', num_lines_overlap)

num_lines_bed = 0

# determine how many lines in original bed file
for line in bed:
    num_lines_bed +=1
print ('Number of total reads = ', num_lines_bed)

# print the fraction of peaks in aligned reads
print ('Fraction of Reads in Peak = ', str(float(num_lines_overlap)/float(num_lines_bed)))
