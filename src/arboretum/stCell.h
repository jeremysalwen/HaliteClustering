#ifndef __STCELL_H
#define __STCELL_H

#include "stCellId.h"

class stCell {
   public:
      stCell() {
	     usedCell = 0;
		 sumOfPoints = 0;
		 id = 0;
		 for (int i=0; i<DIM; i++) {
		    P[i] = 0;
		 }
      }
      ~stCell() { 
		 if (id) {
		    delete id;
		 }		  
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
	  void setId(stCellId *id) {
		 this->stCell::~stCell(); //clean any previous id
		 this->id = new stCellId(); //create a new id
	     *(this->id) = *id; //copy content
	  }
	  void insertPointPartial(stCellId *sonsCellId) {
 	     for (int i=0; i<DIM; i++) {
	        if (!sonsCellId->getBitValue(i)) {
		       P[i]++; // counts the point, since it is in 
			           // this cell's lower half regarding i
	        }
	     }
	  }
	  void copy(stCell *cell) { 
		  //copy the content of this to cell
	  	  cell->usedCell = usedCell;
		  for (int i=0; i<DIM; i++) {
			  cell->P[i] = P[i];
		  } 
		  cell->sumOfPoints = sumOfPoints;
		  cell->setId(id);
	  }

   private:
      char usedCell;
      int P[DIM];      
      int sumOfPoints;
	  stCellId *id;
};

#endif //__STCELL_H
