#!/bin/bash
./Calculate.exe -etas 0.16 -Q 0 -yQ 2 -Nch 1940 -NSamples 50120000 -area 104 > LHCb_eta016_y=2_Q0.txt
./Calculate.exe -etas 0.16 -Q 1 -yQ 2 -Nch 1940 -NSamples 50120000 -area 104 > LHCb_eta016_y=2_Q1.txt
./Calculate.exe -etas 0.32 -Q 0 -yQ 2 -Nch 1940 -NSamples 50120000 -area 104 > LHCb_eta032_y=2_Q0.txt
./Calculate.exe -etas 0.32 -Q 1 -yQ 2 -Nch 1940 -NSamples 50120000 -area 104 > LHCb_eta032_y=2_Q1.txt
