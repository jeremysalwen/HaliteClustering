#ifndef __POINTSOURCE_H
#define __POINTSOURCE_H

/*
* 	Rather than create some interface for containers and iterators, etc, we just have this class
* which is an iterator through the data points, but also has information about the end of the list and the beginning
*/

#include <fstream>
#include <sstream>

namespace Halite {
  template <typename T>
  class PointSource {
  public:
    virtual size_t dimension()      = 0;
    virtual void restartIteration() = 0;
    virtual bool hasNext()   = 0;
    virtual const T* readPoint() = 0;
    virtual ~PointSource() {};
  };
  template <typename T>
    class PackedArrayPointSource : public PointSource<T> {
  public:
    PackedArrayPointSource(const T* array, size_t dim, size_t size) {
      this->data=array;
      this->dim=dim;
      this->size=size;
      this->index=0;
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
    const T* readPoint() {
      return &data[index++*dim];
    }
  private:
    const T* data;
    size_t index;
    size_t dim;
    size_t size;
  };

template <typename T>
  class ArrayOfPointersPointSource : public PointSource<T> {
  public:
    ArrayOfPointersPointSource(T** array, size_t dim, size_t size) {
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
    const T* readPoint() {
      return data[index++];
    }
  private:
    T** data;
    size_t index;
    size_t dim;
    size_t size;
  };
template <typename T>
  class TextFilePointSource : public PointSource<T> {
  public:
  TextFilePointSource(const char* filename):database(filename) {
      if(!database.is_open()) {
	throw std::exception();
      }
      std::getline(database,nextline);
			
      std::istringstream ss(nextline);
      dim=0;
      while(true) {
	T tmp;
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
    const T* readPoint() {
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
    std::vector<T> tmparray;
  };
}
#endif //__POINTSOURCE_H
