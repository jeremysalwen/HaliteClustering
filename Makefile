CXX := g++
SRC := $(wildcard src/*.cpp) $(wildcard src/arboretum/*.cpp) 
OBJ := $(SRC:.cpp=.o)

CXXFLAGS := `pkg-config --cflags opencv` -Wall -Wextra -std=c++11 -g -pedantic
LDFLAGS := -ldb_cxx -lboost_system -lboost_filesystem  `pkg-config --libs opencv` -g

Halite: $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o Halite

demo:   Halite
	./Halite 1e-10 4 1 1 0

clean:
	rm -f ./results/result12d.dat

spotless: clean
	rm -f Halite

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
