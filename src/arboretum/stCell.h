#ifndef __STCELL_H
#define __STCELL_H

#include "stCellId.h"

class stCell {
   public:
     	static stCell* create(int DIM)
    	{
          char* buf = new char[sizeOf(DIM)];
          return ::new (buf) stCell(DIM); 
    	}
        ~stCell() { 
	  if (id) {
	    delete id;
	  }		  
        }
	void destroy()
    	{
		this->~stCell();
       		delete[] (char*)this;
    	}
	static size_t sizeOf(int DIM) {
		return sizeof(stCell) + (DIM - 1) * sizeof(int);
	}
	  void insertPoint() {
	     sumOfPoints++;		 
	  }
	  int getSumOfPoints() {
	     return sumOfPoints;
	  }
	  char getUsedCell() {
	     return usedCell;
	  }	  
	  void useCell() {
	     usedCell = 1;
	  }
	  int getP(int i) {
	     return P[i];
	  }
	  stCellId *getId() {
	     return id;
	  }
	  void setId(stCellId *id,int DIM) {
		 this->stCell::~stCell(); //clean any previous id
		if(id) {
		    this->id = new stCellId(DIM); //create a new id
	     	    *(this->id) = *id; //copy content
		} else {
			this->id=NULL;
		}
	  }
	  void insertPointPartial(stCellId *sonsCellId,int DIM) {
 	     for (int i=0; i<DIM; i++) {
	        if (!sonsCellId->getBitValue(i,DIM)) {
		       P[i]++; // counts the point, since it is in 
			           // this cell's lower half regarding i
	        }
	     }
	  }
	  void copy(stCell *cell,int DIM) { 
		  //copy the content of this to cell
	  	  cell->usedCell = usedCell;
		  for (int i=0; i<DIM; i++) {
			  cell->P[i] = P[i];
		  } 
		  cell->sumOfPoints = sumOfPoints;
		  cell->setId(id,DIM);
	  }

   private:
 	stCell(int DIM) {
	     	usedCell = 0;
		sumOfPoints = 0;
		id = 0;
		for (int i=0; i<DIM; i++) {
                   P[i] = 0;
        	}
     	}

      char usedCell;     
      int sumOfPoints;
      stCellId *id;
      int P[1];
};

#endif //__STCELL_H
