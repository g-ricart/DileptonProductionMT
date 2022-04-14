#!/bin/bash
./Calculate.exe -etas 0.16 -Q 0 -yQ 0 -Nch 2240 -NSamples 50120000 -area 104 > ALICE_eta016_y=0_Q0.txt
./Calculate.exe -etas 0.16 -Q 1 -yQ 0 -Nch 2240 -NSamples 50120000 -area 104 > ALICE_eta016_y=0_Q1.txt
./Calculate.exe -etas 0.32 -Q 0 -yQ 0 -Nch 2240 -NSamples 50120000 -area 104 > ALICE_eta032_y=0_Q0.txt
./Calculate.exe -etas 0.32 -Q 1 -yQ 0 -Nch 2240 -NSamples 50120000 -area 104 > ALICE_eta032_y=0_Q1.txt
