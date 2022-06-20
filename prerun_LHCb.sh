#!/bin/bash

etas=(2.125 2.375 2.625 2.875 3.125 3.375 3.625 3.875 4.125 4.375 4.625 4.875)
multiplicities=(1747 1678 1637 1577 1536 1509 1433 1351 1256 1213 1176 1142)

for (( i = 0; i < 12; i++ )); do
    eta="${etas[$i]}"
    mult="${multiplicities[$i]}"
    file="output/LHCb_$eta.txt"

    echo "./Calculate_dqT.exe -etas 0.32 -Q 1 -yQ 0 -Nch $mult -NSamples 50120000 -area 104 > $file"
done
