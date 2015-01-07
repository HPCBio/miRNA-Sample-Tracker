#!/bin/bash

echo -e "Sample\tRaw reads\tKept Reads\tDropped Reads"
grep 'Input Reads' *.summary.log | perl -pe 's/^(.*)\.trimmomatic\.summary\.log:Input Reads: (\d+) Surviving: (\d+) \(\d+\.\d+%\) Dropped: (\d+) \(\d+\.\d+%\)/$1\t$2\t$3\t$4/'
