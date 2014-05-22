Halite: src/*.cpp src/*.h src/arboretum/*.h
	g++ src/Halite.cpp src/Utile.cpp -o Halite -ldb_cxx `pkg-config --cflags --libs opencv` -Wall -Wextra

demo:   Halite
	./Halite 1e-10 4 1 1 0

clean:
	\rm -f ./results/result12d.dat

spotless: clean
	\rm -f Halite
