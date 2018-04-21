.DEFAULT_GOAL := all

all:
	g++ -O3 -fopenmp smallpt.cpp -o smallpt -s
	
clean:
	rm -rf smallpt.exe
	