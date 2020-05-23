#!/bin/bash

head -n 1 jointGwasMc_HDL.txt >> good_rows.txt
awk -F'\t' '($9<9.9e-309 || $9<2.2250738585072014e-308 && $9>0)' < jointGwasMc_HDL.txt >> 
good_rows.txt
