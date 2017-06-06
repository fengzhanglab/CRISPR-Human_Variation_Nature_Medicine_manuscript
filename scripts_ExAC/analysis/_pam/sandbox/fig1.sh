#!/bin/bash

qsub -q short -P varicas -o fig1.o -e fig1.e -b y -cwd -l h_vmem=6g Rscript fig1.R