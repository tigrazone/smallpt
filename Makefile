.DEFAULT_GOAL := all

all:
	g++ -O3 -fopenmp smallpt.cpp -o smallpt -s
	g++ -O3 -fopenmp smallpt-grid.cpp -o smallpt-grid -s
	
clean:
	rm -rf smallpt.exe smallpt-grid.exe
	