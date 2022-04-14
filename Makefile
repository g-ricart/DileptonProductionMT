all:
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main.cpp -o Calculate.exe -O3 -lgsl -lgslcblas
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main_dqT.cpp -o Calculate_dqT.exe -O3 -lgsl -lgslcblas
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main_dqT_Bolt.cpp -o Calculate_dqT_Bolt.exe -O3 -lgsl -lgslcblas
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main_MT.cpp -o Calculate_mT.exe -O3 -lgsl -lgslcblas
