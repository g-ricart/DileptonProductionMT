all:
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main_dqT.cpp -o Calculate_dqT.exe -O3 -lgsl -lgslcblas
