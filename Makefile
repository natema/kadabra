build:
	g++ -fopenmp -std=c++11 -Ofast -Wall -pedantic -g -Iinclude main.cpp src/* -lm -o kadabra
