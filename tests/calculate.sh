#!/bin/bash
cd $1

# - .vcf.gz files have no stochastic content except a header line with time information
# - .tbi files are generically named, no stochastic content, may be md5sum-checked
# - .fin files are generically named, no stochastic content, may be md5sum-checked

# Therefore:
# - Check md5sums for all types of files, sort

echo ".vcf files:"
for v in *.vcf.gz;do gunzip -c $v | grep -v GATKCommandLine | md5sum;done | sort -V

echo ".tbi files:"
find . -name "*.tbi" | xargs md5sum | sort

echo ".fin files:"
find . -name "*.fin" | xargs md5sum | sort

