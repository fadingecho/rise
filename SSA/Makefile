PARA = -std=c++11 -Wall -O3 
com: rwgraph.cpp SSA_copy.cpp el2bin.cpp option.cpp GLib.hpp
	g++ -c GLib.hpp -o GLib.o $(PARA)
	g++ -c mappedHeap.hpp -o mappedHeap.o $(PARA)
	g++ -c HeapData.hpp -o HeadData.o $(PARA)
	g++ -c option.cpp -o option.o $(PARA)
	g++ -c rwgraph.cpp -o rwgraph.o $(PARA)
	g++ SSA_copy.cpp rwgraph.o option.o -o SSA_copy -fopenmp $(PARA) sfmt/SFMT.c
	g++ el2bin.cpp -o el2bin $(PARA)