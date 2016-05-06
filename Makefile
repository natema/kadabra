build:
	g++ -fopenmp -std=c++11 -Ofast -Wall -g -Iinclude main.cpp src/* -lm -o kadabra
