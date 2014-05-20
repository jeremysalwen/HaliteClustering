MrCC: 
	env DYLD_LIBRARY_PATH="/usr/local/BerkeleyDB.5.0/lib" g++ MrCC.cpp -o MrCC -ldb_cxx-5.0 -I/usr/local/BerkeleyDB.5.0/include/ -L/usr/local/BerkeleyDB.5.0/lib/ `pkg-config --cflags --libs opencv`

demo:   MrCC
	./MrCC 1e-10 4 1 1

clean:
	\rm -f ./results/result12d.dat

spotless: clean
	  \rm -f MrCC
