CXX := g++
SRC := $(wildcard src/*.cpp) $(wildcard src/arboretum/*.cpp) 
OBJ := $(SRC:.cpp=.o)

CXXFLAGS := -Iinclude `pkg-config --cflags opencv` -O2 -g -Wall -Wextra -std=c++11 -pedantic
LDFLAGS := -ldb_cxx -lboost_system -lboost_filesystem  `pkg-config --libs opencv` -O2 -g

Halite: $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o Halite

demo:   Halite
	./Halite 1e-10 4 1 1 0 databases/12d.dat results/result12.dat

clean:
	rm $(OBJ)

spotless: clean
	rm -f Halite

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
