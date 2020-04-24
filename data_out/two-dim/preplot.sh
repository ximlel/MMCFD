#!/bin/bash  --login
shopt -s expand_aliases
alias preplot='~/Softwares/tecplot360ex/bin/preplot'
for i in $(ls *.tec)
do
preplot $i
done
rm *.tec
