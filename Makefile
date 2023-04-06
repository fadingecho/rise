# https://www.jianshu.com/p/442e71755643
# https://www.cnblogs.com/dylancao/p/11947107.html
# https://www.cnblogs.com/jt2001/p/5198733.html
# https://blog.csdn.net/liang_baikai/article/details/110137374

CXX := g++-10
CC := g++-10

WARN_FLAGS := -Wall -Wextra -pedantic
DEBUG_FLAGS := -std=c++20 -Og -ggdb $(WARN_FLAGS)
FLAGS := -std=c++20 -Wall -O3

moga_test: SSA/rwgraph.cpp SSA/el2bin.cpp SSA/option.cpp SSA/GLib.hpp
	$(CXX) -c SSA/GLib.hpp -o GLib.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/mappedHeap.hpp -o mappedHeap.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/HeapData.hpp -o HeadData.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/option.cpp -o option.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/rwgraph.cpp -o rwgraph.o $(DEBUG_FLAGS)
	$(CXX) $(DEBUG_FLAGS) ga_im.cpp rwgraph.o option.o -o moga_demo SSA/sfmt/SFMT.c

soga_test: SSA/rwgraph.cpp SSA/el2bin.cpp SSA/option.cpp SSA/GLib.hpp
	$(CXX) -c SSA/GLib.hpp -o GLib.o $(FLAGS)
	$(CXX) -c SSA/mappedHeap.hpp -o mappedHeap.o $(FLAGS)
	$(CXX) -c SSA/HeapData.hpp -o HeadData.o $(FLAGS)
	$(CXX) -c SSA/option.cpp -o option.o $(FLAGS)
	$(CXX) -c SSA/rwgraph.cpp -o rwgraph.o $(FLAGS)
	$(CXX) $(FLAGS) int_encoded_soga.cpp rwgraph.o option.o -o soga_demo SSA/sfmt/SFMT.c

moga_dbg: SSA/rwgraph.cpp SSA/el2bin.cpp SSA/option.cpp SSA/GLib.hpp
	$(CXX) -c SSA/GLib.hpp -o GLib.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/mappedHeap.hpp -o mappedHeap.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/HeapData.hpp -o HeadData.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/option.cpp -o option.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/rwgraph.cpp -o rwgraph.o $(DEBUG_FLAGS)
	$(CXX) $(DEBUG_FLAGS) int_encoded_copy.cpp rwgraph.o option.o -o moga_debug SSA/sfmt/SFMT.c

moga2: SSA/rwgraph.cpp SSA/el2bin.cpp SSA/option.cpp SSA/GLib.hpp
	$(CXX) -c SSA/GLib.hpp -o GLib.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/mappedHeap.hpp -o mappedHeap.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/HeapData.hpp -o HeadData.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/option.cpp -o option.o $(DEBUG_FLAGS)
	$(CXX) -c SSA/rwgraph.cpp -o rwgraph.o $(DEBUG_FLAGS)
	$(CXX) $(DEBUG_FLAGS) moga.cpp rwgraph.o option.o -o moga2 SSA/sfmt/SFMT.c

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) SSA/*.o