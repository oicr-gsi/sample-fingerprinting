#!/bin/bash
cd $1

find .  -name "*\.vcf.gz" -exec unzip -q {} \; >/dev/null # unzip the results files

# - .vcf files have no stochastic content except a header line with time information
# - .tbi files are generically named, no stochastic content, may be md5sum-checked
# - .fin files are generically named, no stochastic content, may be md5sum-checked

# Therefore:
# - Check md5sums for all types of files, sort

echo ".vcf files:"
find . -name "*.vcf" | xargs grep -v GATKCommandLine | md5sum | sort

echo ".tbi files:"
find . -name "*.tbi" | xargs md5sum | sort

echo ".fin files:"
find . -name "*.fin" | xargs md5sum | sort

