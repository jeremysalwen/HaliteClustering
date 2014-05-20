/**********************************************************************
* GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
*
*                           Homepage: http://gbdi.icmc.usp.br/arboretum
**********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * Ciências Matemáticas e de Computação, University of São Paulo -
 * Brazil (the Databases and Image Group - Intitute of Matematical and
 * Computer Sciences).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *       "This product includes software developed by Grupo de Bases
 *        de Dados e Imagens, Instituto de Ciências Matemáticas e de
 *        Computação, University of São Paulo - Brazil (the Databases 
 *        and Image Group - Intitute of Matematical and Computer
 *        Sciences)"
 *
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names of the research group, institute, university, authors
 *    and collaborators must not be used to endorse or promote products
 *    derived from this software without prior written permission.
 *
 * 5. The names of products derived from this software may not contain
 *    the name of research group, institute or university, without prior
 *    written permission of the authors of this software.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OF THIS SOFTWARE OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * ====================================================================
 *                                            http://gbdi.icmc.usp.br/
 */
 /**
 * @file
 * This file defines the class stCellId.
 *
 * @version 3.0
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
 * @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * 
 */
 // Copyright (c) 2002-2009 GBDI-ICMC-USP

#ifndef __STCELLID_H
#define __STCELLID_H

#include "stCommon.h"
#include <memory.h>
#include <math.h>

//----------------------------------------------------------------------------
// class stCellId
//----------------------------------------------------------------------------
/**
* To identify cells in datasets.
*
* @version 3.0
* @author Caetano Traina Jr (caetano@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
* @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
*/
class stCellId {   

   private:

      /**
      * Vector that is going to be the index in which each bit refers to a dimension.
	  * Possible unused bits will be in the "left side" of index[0].
      */
      unsigned char *index; 

      /**
      * One position for each 8 dimensions.	  
      */
      int nPos;

   public:

      /**
      * Constructor
      */
      stCellId(){		 
         nPos = (int) ceil((double)DIM/8);
         index = new unsigned char[nPos];  // one position for each 8 dimensions
         memset(index,0,nPos); // clean the used memory
      }//end stCellId

      /**
      * Destructor
      */
      ~stCellId(){         
         delete[] index;
      }//end ~stCellId

      /**
      * Gets the value of the bit in position i
      *
      * @param i Bit position.
	  * @return the value of the bit in position i.
      */
      char getBitValue(int i) {		
		int p; // char that stores the bit i
		i += (nPos*8)-DIM; // i + number of unused bits
	    for (p=0; i>7; i-=8,p++);	
		return (index[p] << (i+24) >> 31);
      }//end getBitValue

      /**
      * Inverts the value of the bit in position i
      *
      * @param i Bit position.
      */
      void invertBit(int i) {	    
		int p; // char that stores the bit i
		i += (nPos*8)-DIM; // i + number of unused bits
	    for (p=0; i>7; i-=8,p++);		
        if (index[p] << (i+24) >> 31) {		
		  index[p] -= ( 1 << (7-i) ); // set bit i to 0
        } else {
          index[p] += ( 1 << (7-i) ); // set bit i to 1
        }//endif
      }//end invertBit

      /**
      * Operator = assign a value to a variable.
      *
      * @param cell The stCellId to be assigned.
      */
      void operator = (stCellId &cell) {         
         memcpy(this->index, (static_cast<stCellId &>(cell)).index, nPos);
      }//end operator =

      /**
      * Operator == Compare two values.
      *
      * @param cell The stCellId to be compared.
      * @return 0 if equal, <0 if this < cell, >0  if this > cell.
      */
      int operator == (stCellId &cell) {         
         return (memcmp (this->index, (static_cast<stCellId &>(cell)).index, nPos));
	  }//end operator ==
	
	  void getIndex(unsigned char *index) {
		 memcpy(index, this->index, nPos);
	  }
	
	  void setIndex(unsigned char *index) {
		 memcpy(this->index, index, nPos);
	  }	

};

#endif //__STCELLID_H
