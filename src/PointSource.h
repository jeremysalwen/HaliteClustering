#ifndef __POINTSOURCE_H
#define __POINTSOURCE_H

/*
* 	Rather than create some interface for containers and iterators, etc, we just have this class
* which is an iterator through the data points, but also has information about the end of the list and the beginning
*/
class PointSource {
	public:
		virtual size_t dimension()=0;
		virtual void restartIteration() =0;
		virtual bool endOfIteration() = 0;
		virtual void readPoint(double*) =0;
};

class ArrayOfPointersPointSource : PointSource {
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
		bool endOfIteration() {
			return index>=size;
		}
		void readPoint(double* point) {
			double* thePoint=data[index];
			for(size_t d=0; d<dim; d++) {
				point[d]=thePoint[d];
			}
			index++;
		}
	private:
		double** data;
		size_t index;
		size_t dim;
		size_t size;
};

class TextFilePointSource : PointSource {
	public:
		TextFilePointSource(FILE* file) {
			this->database=file;
			char* line=0;
			size_t n=0;
   			getline(&line, &n, file);
			
			dim=0;
			int errcode;
			do {
				double tmp;
				errcode=fscanf(database,"%lf",&tmp);
				dim++;
			} while(errcode!=EOF);

			dim--; //Last line we do not count because it's the labelling
			free(line);
		}
		size_t dimension()  {
			return dim;
		}
		void restartIteration() {
			fseek(database,0,SEEK_SET);
		}
		bool endOfIteration() {
			int c;
			c=fgetc(database);
			ungetc(c,database);
			return c == EOF;
		}
		void readPoint(double* point) {
			int id;
			for (int j=0; j<dim; j++) {
				fscanf(database,"%lf",&point[j]);
			}//end for
			fscanf(database,"%d",&id); // discarts the class id
		}
	private:
		FILE* database;
		size_t dim;
};

#endif //__POINTSOURCE_H
