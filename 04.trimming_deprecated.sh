#!/bin/bash

#This script is USELESS IF you used cutadapt
for file in *;
       do
                fastx_trimmer -Q33 -l 85 -i $file -o $file.trim
                echo "fastxtrimmer : $file: ok";
        done

for file in *.;
        do
                fastx_trimmer -Q33 -l 85 -i $file -o $file.trim
                echo "fastxtrimmer : $file: ok";
        done

