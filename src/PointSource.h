#ifndef __POINTSOURCE_H
#define __POINTSOURCE_H

/*
* 	Rather than create some interface for containers and iterators, etc, we just have this class
* which is an iterator through the data points, but also has information about the end of the list and the beginning
*/

#include <fstream>
#include <sstream>

class PointSource {
	public:
		virtual size_t dimension()      = 0;
		virtual void restartIteration() = 0;
		virtual bool hasNext()   = 0;
		virtual const double* readPoint() = 0;
		virtual ~PointSource() {};
};
class PackedArrayPointSource : public PointSource {
	public:
		PackedArrayPointSource(double* array, size_t dim, size_t size) {
			this->data=array;
			this->dim=dim;
			this->size=size;
		}
		size_t dimension()  {
			return dim;
		}
		void restartIteration() {
			index=0;
		}
		bool hasNext() {
			return index<size;
		}
		const double* readPoint() {
			return &data[index*dim];
		}
	private:
		double* data;
		size_t index;
		size_t dim;
		size_t size;
};


class ArrayOfPointersPointSource : public PointSource {
	public:
		ArrayOfPointersPointSource(double** array, size_t dim, size_t size) {
			this->data=array;
			this->dim=dim;
			this->size=size;
		}
		size_t dimension()  {
			return dim;
		}
		void restartIteration() {
			index=0;
		}
		bool hasNext() {
			return index<size;
		}
		const double* readPoint() {
			return data[index++];
		}
	private:
		double** data;
		size_t index;
		size_t dim;
		size_t size;
};

class TextFilePointSource : public PointSource {
	public:
		TextFilePointSource(const char* filename):database(filename) {
			if(!database.is_open()) {
				throw std::exception();
			}
   			std::getline(database,nextline);
			
			std::istringstream ss(nextline);
			dim=0;
			while(true) {
				double tmp;
				ss >> tmp;
				if(!ss) break;
				dim++;
			}

			dim--; //Last line we do not count because it's the labelling
			tmparray.resize(dim);
		}
		size_t dimension()  {
			return dim;
		}
		void restartIteration() {
			database.clear();
			database.seekg(0,std::ios_base::beg);
			std::getline(database,nextline);
		}
		bool hasNext() {
			return !database.eof();
		}
		const double* readPoint() {
			std::istringstream ss(nextline);

			for (size_t j=0; j<dim; j++) {
				ss >> tmparray[j];
			}//end for
			int id;
			ss>>id;

			std::getline(database,nextline);
			return &tmparray[0];
		}
	private:
		std::string nextline;
		std::ifstream database;
		size_t dim;
		std::vector<double> tmparray;
};

#endif //__POINTSOURCE_H
