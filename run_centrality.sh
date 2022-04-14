#!/bin/bash
./Calculate.exe -etas 0.32 -Q 1 -Nch 1770 -area 96 -NSamples 50120000 > 0-10_eta032_y=2_Q1.txt
./Calculate.exe -etas 0.32 -Q 1 -Nch 1181 -area 71 -NSamples 50120000 > 10-20_eta032_y=2_Q1.txt
./Calculate.exe -etas 0.32 -Q 1 -Nch 792 -area 54 -NSamples 50120000 > 20-30_eta032_y=2_Q1.txt
./Calculate.exe -etas 0.32 -Q 1 -Nch 514 -area 41 -NSamples 50120000 > 30-40_eta032_y=2_Q1.txt

