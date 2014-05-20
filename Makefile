MrCC: 
	g++ src/Halite.cpp -o Halite -ldb_cxx `pkg-config --cflags --libs opencv`

demo:   Halite
	./Halite 1e-10 4 1 1 12 0

clean:
	\rm -f ./results/result12d.dat

spotless: clean
	\rm -f Halite
