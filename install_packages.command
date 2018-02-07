#!/bin/bash

cd "${0%/*}"
cd src

Rscript installpackages.R

echo
echo
echo
echo \>\>\> Press Enter to finish
#read pressed