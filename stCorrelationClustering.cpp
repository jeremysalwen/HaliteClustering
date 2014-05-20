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
* This file implements the class stCorrelationClustering.
*
* @version 1.0
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Caetano Traina Jr (caetano@icmc.usp.br)
*/
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// class stCorrelationClustering
//------------------------------------------------------------------------------

#include "stCorrelationClustering.h"

stCorrelationClustering::stCorrelationClustering (double **objectsArray, FILE *database, int normalizeFactor, int centralConvolutionValue, 
												  int neighbourhoodConvolutionValue, double pThreshold, int H, int hardClustering, 
												  int initialLevel, DBTYPE dbType, char memory) {

    // stores H, hardClustering and initialLevel
    this->H = H;
	this->hardClustering = hardClustering;
	this->initialLevel = initialLevel;

    // builds the counting tree and inserts objects on it
    calcTree = new stCountingTree(H, dbType, (memory==2));
    fastDistExponent(objectsArray, database, normalizeFactor, memory);

    // builds vectors to describe the positions of a cluster center in the data space
    minBetaClusterCenter = new double[DIM];
    maxBetaClusterCenter = new double[DIM];

    // builds vectors to the parents of a cluster center and of a neighbour
    betaClusterCenterParents = new stCell*[H];
    neighbourParents = new stCell*[H];
	for (int i=0; i<H; i++) {
		betaClusterCenterParents[i] = new stCell();
        neighbourParents[i] = new stCell();
	}

    // builds auxiliar vectors used to search for the relevant attributes
    attributesRelevance = new double[DIM];

    // builds auxiliar vector to describe which neighbours belong to a found cluster
    neighbourhood = new char[DIM];

    // stores the convolution matrix (center and direct neighbours)
    this->centralConvolutionValue=centralConvolutionValue;
    this->neighbourhoodConvolutionValue=neighbourhoodConvolutionValue;

    // stores the pThreshold
    this->pThreshold = pThreshold;

    // initiates the number of found clusters
    numBetaClusters=numCorrelationClusters=0;

    // defines the maximum number of beta-clusters
    int maxNumBetaClusters=2000;

    // builds pointers to describe the found clusters
    minBetaClusters = new double*[maxNumBetaClusters];
    maxBetaClusters = new double*[maxNumBetaClusters];
    dimBetaClusters = new char*[maxNumBetaClusters];
	dimCorrelationClusters = new char*[maxNumBetaClusters];
	correlationClustersBelongings = new int[maxNumBetaClusters];
	costBetaClusters = new int[maxNumBetaClusters];	
	levelBetaClusters = new int[maxNumBetaClusters];		
    for (int i=0; i<maxNumBetaClusters; i++) { // initiation
      minBetaClusters[i]=maxBetaClusters[i]=0;
      dimBetaClusters[i]=dimCorrelationClusters[i]=0;
	  correlationClustersBelongings[i]=costBetaClusters[i]=levelBetaClusters[i]=-1;
    }//end for
	
	//create memory space for neighbour and betaClusterCenter
	neighbour = new stCell();
	betaClusterCenter = new stCell();

}//end stCorrelationClustering::stCorrelationClustering

//---------------------------------------------------------------------------
stCorrelationClustering::~stCorrelationClustering() {

    // disposes the used structures
	delete neighbour;
	delete betaClusterCenter;
    delete [] neighbourhood;    
    delete [] attributesRelevance;
    delete [] minBetaClusterCenter;
    delete [] maxBetaClusterCenter;
	delete [] correlationClustersBelongings;
	for (int i=0; i<H; i++) {
		delete betaClusterCenterParents[i];
        delete neighbourParents[i];
	}
    delete [] betaClusterCenterParents;
    delete [] neighbourParents;
	for (int i=0; i<numBetaClusters; i++) {
		delete [] minBetaClusters[i];
		delete [] maxBetaClusters[i];
		delete [] dimBetaClusters[i];
	}//end for
	for (int i=0; i<numCorrelationClusters; i++) {
		delete [] dimCorrelationClusters[i];	
	}//end for
	delete [] dimCorrelationClusters;
    delete [] minBetaClusters;
    delete [] maxBetaClusters;
    delete [] dimBetaClusters;
	delete [] costBetaClusters;
	delete [] levelBetaClusters;	
    delete calcTree;

}//end stCorrelationClustering::~stCorrelationClustering

//---------------------------------------------------------------------------
void stCorrelationClustering::findCorrelationClusters() {

    // defines when a new cluster is found
    int ok, center, total;
    do { // looks for a cluster in each loop
      ok=0; // no new cluster was found
      // defines the initial grid level to analyze
      int level=initialLevel;
      do { // analyzes each level until a new cluster is found
        // apply the convolution matrix to each grid cell in the current level
        // finds the cell with the bigest convolution value		
		if (walkThroughConvolution(level)) {          
          betaClusterCenter->useCell(); // visited cell
		  calcTree->commitCell(betaClusterCenterParents, betaClusterCenter, level); // commit changes in the tree
		  if (level) {
			// pointer to a neighbour of the father
			stCell *fatherNeighbour = new stCell();	
			stCell **cellP = new stCell*[level-1];
			for (int j=0;j<level-1;j++) {
				cellP[j] = new stCell();
			}//end for
			for (int i=0; i<DIM; i++) {			   
	  			// initiates total with the number of points in the father
				total=betaClusterCenterParents[level-1]->getSumOfPoints();
				// discovers the number of points in the center
				if (betaClusterCenter->getId()->getBitValue(i)) {
  				  center = total - betaClusterCenterParents[level-1]->getP(i);
				} else {
				  center = betaClusterCenterParents[level-1]->getP(i);
				}//end if
				
				// looks for the points in the direct neighbours of the father
				if (internalNeighbour(i,betaClusterCenterParents[level-1],&fatherNeighbour,betaClusterCenterParents,level-1)) {
				  total += fatherNeighbour->getSumOfPoints();
				}//end if
				
				//copy parents
				for (int j=0;j<level-1;j++) {
					betaClusterCenterParents[j]->copy(cellP[j]);
				}//end for
				
                // looks for the external neighbour
				if (externalNeighbour(i,betaClusterCenterParents[level-1],&fatherNeighbour,betaClusterCenterParents,cellP,level-1)) {
				  total += fatherNeighbour->getSumOfPoints();
				}//end if
				
  			    // percentual of points in the center related to the average
                attributesRelevance[i] = (100*center)/((double)total/6);				
				// right critical value for the statistical test
				int criticalValue = GetCriticalValueBinomialRight2(total, (double)1/6, pThreshold);
				if (center > criticalValue) {
				  ok=1; // new cluster found
				}//end if
			}//end for
			delete fatherNeighbour;
			for (int j=0;j<level-1;j++) {
				delete cellP[j];
			}//end for
			delete [] cellP;
		  } else { // analyzes each dimension based on the points distribution of the entire database
				// initiate the total of points
				total=calcTree->getSumOfPoints();			
				for (int i=0; i<DIM; i++) {			   
					// discovers the number of points in the center
					if (betaClusterCenter->getId()->getBitValue(i)) {
						center = total - calcTree->getP()[i];
					} else {
						center = calcTree->getP()[i];
    				}//end if			
  			        // percentual of points in the center related to the average
                    attributesRelevance[i] = (100*center)/((double)total/2);					
					// right critical value for the statistical test
					int criticalValue = GetCriticalValueBinomialRight2(total, (double)1/2, pThreshold);
					if (center > criticalValue) {
					  ok=1; // new cluster found
					}//end if
				}//end for
		  }//end if		
		}//end if

        if (!ok) {
          level++; // next level to be analyzed
        }//end if
      } while (!ok && level < H);//end do while

      if (ok) { // if a new cluster was found...	    

        // discovers the cThreshold based on the minimum description length method
        double cThreshold = calcCThreshold(attributesRelevance);

        // new cluster found
        numBetaClusters++;
		printf("a beta-cluster was found at the Counting-tree level %d.\n",level); // prints the level in which a new beta-cluster was found
		levelBetaClusters[numBetaClusters-1] = level;
        minBetaClusters[numBetaClusters-1] = new double[DIM];
        maxBetaClusters[numBetaClusters-1] = new double[DIM];
        dimBetaClusters[numBetaClusters-1] = new char[DIM];

        // important dimensions
        for (int i=0;i<DIM;i++) {
          dimBetaClusters[numBetaClusters-1][i]=(attributesRelevance[i] >= cThreshold);
        }//end for

        // analyzes neighbours in important dimensions to verify which of them also belong to the found cluster
        for (int i=0;i<DIM;i++) {
          neighbourhood[i]='N'; // no direct neighbour belongs to the cluster
        }//end for
        // center's position in the data space
        cellPosition(betaClusterCenter,betaClusterCenterParents,minBetaClusterCenter,maxBetaClusterCenter,level);
        // for each important dimension, analyzes internal and external neighbours to decide if they also
        // belong to the cluster
        for (int i=0;i<DIM;i++) {
          if (dimBetaClusters[numBetaClusters-1][i]) {
            // looks for the internal neighbour
            if (internalNeighbour(i,betaClusterCenter,&neighbour,betaClusterCenterParents,level)) { // internal neighbour in important dimension always belongs to the cluster
              // neighbour's position in the data space
              cellPositionDimensionE_j(neighbour,betaClusterCenterParents,&minNeighbour,&maxNeighbour,level,i);
              if (maxBetaClusterCenter[i] > maxNeighbour) {
                neighbourhood[i]='I'; // inferior neighbour in i belongs to the cluster
              } else {
                neighbourhood[i]='S'; // superior neighbour in i belongs to the cluster
              }//end if
            }//end if

            // looks for the external neighbour
            for (int j=0;j<level;j++) {
			    betaClusterCenterParents[j]->copy(neighbourParents[j]);	
			}//end for
            if (externalNeighbour(i,betaClusterCenter,&neighbour,betaClusterCenterParents,neighbourParents,level)) { // analyzes external neighbour to decide if it belongs to the cluster
              if (neighbourhood[i] == 'N') {
                // neighbour's position in the data space
                cellPositionDimensionE_j(neighbour,neighbourParents,&minNeighbour,&maxNeighbour,level,i);
                if (maxBetaClusterCenter[i] > maxNeighbour) {
                  neighbourhood[i]='I'; // inferior neighbour in i belongs to the cluster
                } else {
                  neighbourhood[i]='S'; // superior neighbour in i belongs to the cluster
                }//end if
              } else {
                neighbourhood[i]='B'; // both inferior and superior neighbours in i belong to the cluster
              }//end if
            }//end if
          }//end if
        }//end for        

        // stores the description of the found cluster
        double length = maxBetaClusterCenter[0]-minBetaClusterCenter[0];
        for (int i=0;i<DIM;i++) {
          if (dimBetaClusters[numBetaClusters-1][i]) { // dimension important to the cluster
            // analyzes if the neighbours in i also belong to the cluster
            switch (neighbourhood[i]) {
              case 'B': // both inferior and superior neighbours in i belong to the cluster
                  minBetaClusterCenter[i]-=length;
                  maxBetaClusterCenter[i]+=length;
                  break;
              case 'S': // superior neighbour in i belongs to the cluster
                  maxBetaClusterCenter[i]+=length;
                  break;
              case 'I': // inferior neighbour in i belongs to the cluster
                  minBetaClusterCenter[i]-=length;
            }//end switch
            // new cluster description - relevant dimension
            minBetaClusters[numBetaClusters-1][i] = minBetaClusterCenter[i];
            maxBetaClusters[numBetaClusters-1][i] = maxBetaClusterCenter[i];
          } else {
            // new cluster description - irrelevant dimension
            minBetaClusters[numBetaClusters-1][i] = 0;
            maxBetaClusters[numBetaClusters-1][i] = 1;
          }//end if
        }//end for
      }//end if
      // stops when no new cluster is found
    } while (ok);//end do while

	printf("\n%d beta-clusters were found.\n",numBetaClusters); // prints the number of beta-clusters found
	mergeBetaClusters(); // merges clusters that share some database space
	printf("\n%d correlation clusters left after the merging fase.\n",numCorrelationClusters); // prints the number of correlation clusters found

}//end stCorrelationClustering::findCorrelationClusters

//---------------------------------------------------------------------------
double stCorrelationClustering::calcCThreshold(double *attributesRelevance) {

  double *sortedRelevance = new double[DIM];
  for (int i=0;i<DIM;i++) {
    sortedRelevance[i]=attributesRelevance[i];
  }//end for  
  qsort(sortedRelevance,DIM,sizeof(double),compare); // sorts the relevances vector
  double cThreshold = sortedRelevance[minimumDescriptionLength(sortedRelevance)];
  delete [] sortedRelevance;
  return cThreshold;

}//end stCorrelationClustering::calcCThreshold

//---------------------------------------------------------------------------
int stCorrelationClustering::minimumDescriptionLength(double *sortedRelevance) {

  int cutPoint=-1;
  double preAverage, postAverage, descriptionLength, minimumDescriptionLength;
  for (int i=0;i<DIM;i++) {
    descriptionLength=0;
    // calculates the average of both sets
    preAverage=0;
    for (int j=0;j<i;j++) {
      preAverage+=sortedRelevance[j];
    }//end for
    if (i) {
      preAverage/=i;      
	  descriptionLength += (ceil(preAverage)) ? (log10(ceil(preAverage))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end if
    postAverage=0;
    for (int j=i;j<DIM;j++) {
      postAverage+=sortedRelevance[j];
    }//end for
    if (DIM-i) {
      postAverage/=(DIM-i);      
	  descriptionLength += (ceil(postAverage)) ? (log10(ceil(postAverage))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end if
    // calculates the description length
    for (int j=0;j<i;j++) {      
	  descriptionLength += (ceil(fabs(preAverage-sortedRelevance[j]))) ? (log10(ceil(fabs(preAverage-sortedRelevance[j])))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end for
    for (int j=i;j<DIM;j++) {      
	  descriptionLength += (ceil(fabs(postAverage-sortedRelevance[j]))) ? (log10(ceil(fabs(postAverage-sortedRelevance[j])))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end for
    // verify if this is the best cut point
    if (cutPoint==-1 || descriptionLength < minimumDescriptionLength) {
      cutPoint=i;
      minimumDescriptionLength = descriptionLength;
    }//end if
  }//end for
  return cutPoint;

}//end stCorrelationClustering::minimumDescriptionLength

//---------------------------------------------------------------------------
int stCorrelationClustering::walkThroughConvolution(int level) {
		
	//try to get the db in the current level
	Db *db = calcTree->getDb(level);
	if (db) { //got it
		// pointers to the parents of a cell
		stCell **parentsVector = new stCell*[level];
		for (int i=0;i<level;i++) {
			parentsVector[i] = new stCell();
		}//end for
	    
		//prepare the fullId array
		int nPos = (int) ceil((double)DIM/8);
		unsigned char *fullId = new unsigned char[(level+1)*nPos];
		unsigned char *ccFullId = new unsigned char[(level+1)*nPos];
		memset(fullId,0,(level+1)*nPos);
		memset(ccFullId,0,(level+1)*nPos);
		
		//prepare the cell and the Dbts to receive 
		//key/data pairs from the dataset 
		stCell *cell = new stCell();
		stCellId *id = new stCellId();
		Dbt searchKey, searchData;
		searchKey.set_data(fullId);
		searchKey.set_ulen((level+1)*nPos);
		searchKey.set_flags(DB_DBT_USERMEM);
		searchData.set_data(cell);
		searchData.set_ulen(sizeof(stCell));
		searchData.set_flags(DB_DBT_USERMEM);
		
		// Get a cursor
		Dbc *cursorp;
		db->cursor(NULL,&cursorp,0);
		
		//prepare to walk through the level and find the
		//cell with the biggest convoluted value
		int biggestConvolutionValue = -MAXINT;
		int newConvolutionValue;
		char clusterFoundBefore;
		double *maxCell = new double[DIM];
		double *minCell = new double[DIM];

		// iterate over the database, retrieving each record in turn
		int ret;
		while ((ret = cursorp->get(&searchKey, &searchData, DB_NEXT)) == 0) {
			// Does not analyze cells analyzed before and cells that can't be the biggest convolution center.
			// It speeds up the algorithm, specially when neighbourhoodConvolutionValue <= 0
			if ( (!cell->getUsedCell()) && ( (neighbourhoodConvolutionValue > 0) || 
											 ((cell->getSumOfPoints()*centralConvolutionValue) > biggestConvolutionValue) ||
											 ( ((cell->getSumOfPoints()*centralConvolutionValue) == biggestConvolutionValue) && 
											   (memcmp(fullId, ccFullId, (level+1)*nPos) < 0) ) ) ) {
				//set id for cell
				id->setIndex(fullId+(level*nPos)); //copy from fullId to id
				cell->setId(id); //copy from id to cell->id
				//finds the parents of cell
				calcTree->findParents(fullId,parentsVector,level);
				// discovers the position of cell in the data space
				cellPosition(cell,parentsVector,minCell,maxCell,level);
				// verifies if this cell belongs to a cluster found before
				clusterFoundBefore=0;
				for (int i=0;!clusterFoundBefore && i<numBetaClusters;i++) {
					clusterFoundBefore = 1;
					for (int j=0; clusterFoundBefore && j<DIM; j++) {
						// Does not cut off cells in a level upper than the level where a cluster was found
						if (!(maxCell[j] <= maxBetaClusters[i][j] && minCell[j] >= minBetaClusters[i][j])) {
							clusterFoundBefore = 0;
						}//end if
					}//end for
				}//end for
				
				if (!clusterFoundBefore) { // the cell doesn't belong to any found cluster
					// applies the convolution matrix to cell
					if (neighbourhoodConvolutionValue) {
						newConvolutionValue=applyConvolution(cell,parentsVector,level);
					} else {
						newConvolutionValue=centralConvolutionValue*(cell->getSumOfPoints()); // when the neighbourhood weight is zero
					}//end if
					if ( (newConvolutionValue > biggestConvolutionValue) || 
						 ((newConvolutionValue == biggestConvolutionValue) && (memcmp(fullId, ccFullId, (level+1)*nPos) < 0))) {
						// until now, cell is the biggest convolution value, thus, set the new biggest value and copy
						// cell and its parents to betaClusterCenter and its parents
						biggestConvolutionValue = newConvolutionValue;
						memcpy(ccFullId, fullId, (level+1)*nPos);
						cell->copy(betaClusterCenter);
						for (int j=0;j<level;j++) {
							parentsVector[j]->copy(betaClusterCenterParents[j]);
						}//end for
					}//end if
				}//end if
			}//end if        
		}//end while
		if (ret != DB_NOTFOUND) { //it should never enter here
			cout << "Error!" << endl;
			return 0; //error
		}
		
		//closes the cursor
		if (cursorp != NULL) {
			cursorp->close();
		}
		
		// disposes the used memory
		delete [] minCell;
		delete [] maxCell;
		delete cell;
		delete id;
		delete [] fullId;
		delete [] ccFullId;
		for (int i=0;i<level;i++) {
			delete parentsVector[i];
		}//end for
		delete [] parentsVector;
		
		return 1; //Success
	}
	return 0; //Error
}//end stCorrelationClustering::walkThroughConvolution

//---------------------------------------------------------------------------
int stCorrelationClustering::applyConvolution(stCell *cell, stCell **cellParents, int level) {

  stCell *neighbour = new stCell();
  stCell **neighbourParents = new stCell*[level];
  for (int j=0;j<level;j++) {
	  neighbourParents[j] = new stCell();
  }//end for
  int newValue = cell->getSumOfPoints()*centralConvolutionValue;
  // looks for the neighbours
  for (int k=0;k<DIM;k++) {
    if (internalNeighbour(k,cell,&neighbour,cellParents,level)  ) {
      newValue+=(neighbour->getSumOfPoints()*neighbourhoodConvolutionValue);
    }//end if
	
	//copy parents
	for (int j=0;j<level;j++) {
	   cellParents[j]->copy(neighbourParents[j]);
	}//end for
	  
    if (externalNeighbour(k,cell,&neighbour,cellParents,neighbourParents,level)) {
		newValue+=(neighbour->getSumOfPoints()*neighbourhoodConvolutionValue);
	}//end if
  }//end for
  delete neighbour;
  for (int j=0;j<level;j++) {
	  delete neighbourParents[j];
  }//end for
  delete [] neighbourParents;	
  // return the cell value after applying the convolution matrix
  return newValue;

}//end stCorrelationClustering::applyConvolution

//---------------------------------------------------------------------------
void stCorrelationClustering::cellPosition(stCell *cell, stCell **cellParents, 
										   double *min, double *max, int level) {

  if (cell) {
    if (level) {
      cellPosition(cellParents[level-1],cellParents,min,max,level-1);
      for (int i=0; i<DIM; i++) {
        if (cell->getId()->getBitValue(i)) { // bit in the position i is 1
          min[i] += ((max[i]-min[i])/2);
        } else { // bit in the position i is 0
          max[i] -= ((max[i]-min[i])/2);
        }//end if
      }//end for
    } else { // level zero
      for (int i=0; i<DIM; i++) {
        if (cell->getId()->getBitValue(i)) { // bit in the position i is 1
          min[i] = 0.5;
          max[i] = 1;
        } else { // bit in the position i is 0
          min[i] = 0;
          max[i] = 0.5;
        }//end if
      }//end for
    }//end if
  }//end if

}//end stCorrelationClustering::cellPosition

//---------------------------------------------------------------------------
void stCorrelationClustering::cellPositionDimensionE_j(stCell *cell, stCell **cellParents, 
													   double *min, double *max, int level, int j) {

  if (cell) {
    if (level) {
      cellPositionDimensionE_j(cellParents[level-1],cellParents,min,max,level-1,j);      
      if (cell->getId()->getBitValue(j)) { // bit in the position j is 1
        *min += ((*max-*min)/2);
      } else { // bit in the position j is 0
        *max -= ((*max-*min)/2);
      }//end if      
    } else { // level zero
      if (cell->getId()->getBitValue(j)) { // bit in the position j is 1
        *min = 0.5;
        *max = 1;
      } else { // bit in the position j is 0
        *min = 0;
        *max = 0.5;
      }//end if
    }//end if
  }//end if

}//end stCorelationClustering::cellPositionDimensionE_j

//---------------------------------------------------------------------------
int stCorrelationClustering::externalNeighbour(int dimIndex, stCell *cell, stCell **neighbour,
											   stCell **cellParents, stCell **neighbourParents, int level) {
   if (level) {
	 int found;
     stCell *father = cellParents[level-1];
     if (cell->getId()->getBitValue(dimIndex) ^ father->getId()->getBitValue(dimIndex)) { // XOR - different bit values -> starts going down in the tree
	   //finds the internal neighbour of the father
	   found = internalNeighbour(dimIndex, father, &neighbourParents[level-1], cellParents, level-1);
     } else {  // equal bit values -> continues going up in the tree
	   // recursively, finds the external neighbour of the father
       found = externalNeighbour(dimIndex,father,&neighbourParents[level-1],cellParents,neighbourParents,level-1); 
     }//end if
     if (found) { // father's neighbour was found
	   // find the external neighbour of cell in dimension i
       return internalNeighbour(dimIndex, cell, neighbour, neighbourParents, level);
     }//end if
     return 0; // there is no father's neighbour
   }//end if
   return 0; // a cell in level zero never has an external neighbour

}//end stCorrelationClustering::externalNeighbour

//---------------------------------------------------------------------------
int stCorrelationClustering::internalNeighbour(int dimIndex, stCell *cell, stCell **neighbour,
											   stCell **cellParents, int level) {

  // creates the id that the neighbour should have
  stCellId *neighboursId = new stCellId();
  *neighboursId=*cell->getId();
  neighboursId->invertBit(dimIndex);
  int found = calcTree->findInNode(cellParents, neighbour, neighboursId, level);
  delete neighboursId;
  return found;

}//end stCorrelationClustering::internalNeighbour

//---------------------------------------------------------------------------
void stCorrelationClustering::fastDistExponent(double **objectsArray, FILE *database, int normalizeFactor, char memory) {

   double *minD, *maxD, biggest;
   double *onePoint, *resultPoint, *a, *b; // y=Ax+B to normalize each dataset.
   double normalizationFactor = 1.0;

   minD = (double *) calloc ((1+DIM),sizeof(double));
   maxD = (double *) calloc ((1+DIM),sizeof(double));
   a = (double *) calloc(DIM,sizeof(double));
   b = (double *) calloc(DIM,sizeof(double));
   onePoint = (double *) calloc(DIM,sizeof(double));
   resultPoint = (double *) calloc(DIM,sizeof(double));

   // normalizes the data
   minMax(objectsArray, database, minD, maxD, memory);
   biggest = 0; // for Normalize==0, 1 or 3
   // normalize=0->Independent, =1->mantain proportion, =2->Clip
   //          =3->Geo Referenced
   if (normalizeFactor == 2) {
     biggest = MAXDOUBLE;
   }//end if

   for (int i=0; i<DIM; i++) {
     a[i] = (maxD[i] - minD[i]) * normalizationFactor; //a takes the range of each dimension
     b[i] = minD[i];
     if (a[i] == 0) {
       a[i] = 1;
     }//end if
   }//end for

   for (int i=0; i<DIM; i++) {
     if ((normalizeFactor < 2 || normalizeFactor == 3) && biggest < a[i]) {
       biggest = a[i];
     }//end if
     if (normalizeFactor == 2 && biggest > a[i]) {
       biggest = a[i];
     }//end if
   }//end for

   if (normalizeFactor != 0) {
     for (int i=0; i<DIM; i++) {
       a[i] = biggest; // normalized keeping proportion
     }//end for
     /* when we have the proportional normalization, every A[i] are gonna receive
     the biggest range.*/
   }//end if

   if (normalizeFactor >= 0) {
      calcTree->setNormalizationVectors(a,b); // if there is some normalization
   }//end if

   //process each point
   if (memory != 0) { // limited or none memory
      fseek(database,0,SEEK_SET); //go to the file begin
   }
   for (int i=0; i<SIZE; i++) {
	 (memory == 0) ? copyPoint(objectsArray[i], onePoint) : readPoint(database, onePoint);
	 calcTree->insertPoint(onePoint,resultPoint); //add to the grid structure
   }//end for

   // disposes used memory
   delete[] onePoint;
   delete[] resultPoint;
   delete[] a;
   delete[] b;
   delete[] minD;
   delete[] maxD;

}//end stCorrelationClustering::FastDistExponent

//---------------------------------------------------------------------------
void stCorrelationClustering::minMax(double **objectsArray, FILE *database, double *min, double *max, char memory) {

  timeNormalization = clock(); //start normalization time
  double *onePoint = new double[DIM];
  for (int j=0; j<DIM; j++){ // sets the values to the minimum/maximum possible here
    min[j] = MAXDOUBLE;
    max[j] = -MAXDOUBLE;
  }// end for
  // looking for the minimum and maximum values
  if (memory != 0) { // limited or none memory
    fseek(database,0,SEEK_SET); //go to the file begin
  }
  for (int i=0; i<SIZE; i++) {
	(memory == 0) ? copyPoint(objectsArray[i], onePoint) : readPoint(database, onePoint);  
    for (int j=0; j<DIM; j++) {
      if (onePoint[j] < min[j]) {
        min[j] = onePoint[j];
      }//end if
      if (onePoint[j] > max[j]) {
        max[j] = onePoint[j];
      }//end if
    }//end for
  }//end for
  delete [] onePoint;
  timeNormalization = (clock()-timeNormalization); //total time spent in the normalization

}//end stCorrelationClustering::MinMax

//---------------------------------------------------------------------------
void stCorrelationClustering::mergeBetaClusters() {

	int i=0, j, k, aux;
	// merges beta-clusters
	while (i<numBetaClusters) {
		j=i+1;		
		while (j<numBetaClusters) {
			if (shouldMerge(i, j)) { // merges both beta-clusters
				if (correlationClustersBelongings[i]==-1 && correlationClustersBelongings[j]==-1) { // both clusters belong to no merged cluster
					correlationClustersBelongings[i]=correlationClustersBelongings[j]=numCorrelationClusters++; // new merged cluster
				} else {
					if (correlationClustersBelongings[i]!=-1 && correlationClustersBelongings[j]!=-1) { // both clusters belong to some merged cluster(s)
						if (correlationClustersBelongings[i]!=correlationClustersBelongings[j]) { // both clusters belong to different merged clusters
							numCorrelationClusters--;
							for (k=0; k<numBetaClusters; k++) {
								if (k!=i && k!=j && (correlationClustersBelongings[k]==correlationClustersBelongings[i] || correlationClustersBelongings[k]==correlationClustersBelongings[j])) {
									correlationClustersBelongings[k]=(correlationClustersBelongings[i]>correlationClustersBelongings[j]) ? correlationClustersBelongings[j] : correlationClustersBelongings[i];
								}//end if
							}//end for
							if (correlationClustersBelongings[i]>correlationClustersBelongings[j]) {
								aux = correlationClustersBelongings[i]; // deleted cluster 
								correlationClustersBelongings[i]=correlationClustersBelongings[j];								 
							} else { 
								aux = correlationClustersBelongings[j]; // deleted cluster 
								correlationClustersBelongings[j]=correlationClustersBelongings[i];
							}//end if
							for (k=0; k<numBetaClusters; k++) {
								if (correlationClustersBelongings[k] > aux) {
									correlationClustersBelongings[k]--;
								}//end if
							}//end for
						}//end if
					} else { // only one of the beta-clusters belongs to some merged cluster
						(correlationClustersBelongings[i]==-1) ? correlationClustersBelongings[i]=correlationClustersBelongings[j] : correlationClustersBelongings[j]=correlationClustersBelongings[i];
					}//end if
				}//end if
			}//end if
			j++; // next beta-cluster
		}//end while
		if (correlationClustersBelongings[i] == -1) {
			correlationClustersBelongings[i] = numCorrelationClusters++; // new merged cluster
		}//end if
		i++; // next cluster
	}//end while

	// important dimensions to the merged clusters
	for (i=0; i<numCorrelationClusters; i++) {
		dimCorrelationClusters[i] = new char[DIM];
		for (j=0; j<DIM;j++) {
			dimCorrelationClusters[i][j]=0;
		}//end for
	}//end for
	for (i=0; i<numBetaClusters; i++) {
		for (j=0; j<DIM; j++) {
			if (!dimCorrelationClusters[correlationClustersBelongings[i]][j]) {
				dimCorrelationClusters[correlationClustersBelongings[i]][j] = dimBetaClusters[i][j];
			}//end if
		}//end for
	}//end for

}//end stCorrelationClustering::mergeBetaClusters

int stCorrelationClustering::shouldMerge(int i, int j) {
	
	// discovers if beta-cluster i shares database space with beta-cluster j
	int shareSpace=1;
	for(int k=0; shareSpace && k<DIM; k++) {
		if (!(maxBetaClusters[i][k] > minBetaClusters[j][k] && minBetaClusters[i][k] < maxBetaClusters[j][k])) {
			shareSpace=0; // beta-clusters i and j do not share database space
		}//end if
	}//end for
	
	if (shareSpace) {
		if (!hardClustering) {
			//does a PCA based analysis
			if (costBetaClusters[i]==-1) {
				costBetaClusters[i] = cost(i,-1);
			}
			if (costBetaClusters[j]==-1) {
				costBetaClusters[j] = cost(j,-1);
			}			
			return ((double)(costBetaClusters[i]+costBetaClusters[j])/cost(i,j)) >= 1; //merges if the merged cluster compacts best
		}
		return 1; //merge
	}
	return 0; //not merge
}

int stCorrelationClustering::cost(int i, int j) {
	
	int clusterSize=0, cost=0;
	
	//prepare the input for PCA
	cv::Mat clusterMat = inputPCA(i, j, &clusterSize);
	
	//applies PCA in the cluster
	if (clusterSize < DIM) {
	    return 0; // not possible to apply PCA
	}
	cv::PCA princomp(clusterMat, // pass the data
					 cv::Mat(), // we do not have a pre-computed mean vector, so let the PCA engine to compute it
					 CV_PCA_DATA_AS_ROW, // indicate that the vectors are stored as matrix rows
					 DIM // specify, how many principal components to retain
					 );
	clusterMat = princomp.project(clusterMat);
	
	//finds the minimum and maximum values of cluster in each PCA axis
	double *min = new double[DIM];
	double *max = new double[DIM];
	for (int d=0; d<DIM; d++) { // sets the values to the minimum/maximum possible here
		min[d] = MAXDOUBLE;
		max[d] = -MAXDOUBLE;
	}// end for
	for (int p=0; p<clusterSize; p++) {
		for (int d=0; d<DIM; d++) {
			if (clusterMat.at<double>(p,d) > max[d]) {
				max[d] = clusterMat.at<double>(p,d);
			}
			if (clusterMat.at<double>(p,d) < min[d]) {
				min[d] = clusterMat.at<double>(p,d);
			}
		}
	}
	
	//cost for the points
	for (int p=0; p<clusterSize; p++) {
		for (int d=0; d<DIM; d++) {
			cost += indCost(clusterMat.at<double>(p,d)-(((max[d]-min[d])/2) + min[d])); //distance to the center in each dimension
		}
	}
	
	//cost for eigenvectors, min and max
	for (int d=0; d<DIM; d++) {
		cost += indCost(min[d]) + indCost(max[d]);
		for (int k=0; k<DIM; k++) {
			cost += indCost(princomp.eigenvectors.at<double>(d,k));
		}
	}	
	
	//disposes the used memory
	clusterMat.release();
	delete[] min;
	delete[] max;
	
	return cost;
}

int stCorrelationClustering::indCost(double n) {
	n = ceil(fabs(n*1000000)); //ignores the sign, since + and - cost the same, and 
	//considers the cost of the smallest integer bigger than n
	if (n <= 1) {
		return 1; // zero and one cost one
	}
	return (int) ceil(log(n)/log((double) 2)); // cost of n, when n > 1
}

cv::Mat stCorrelationClustering::inputPCA(int i, int j, int *clusterSize) {
	//prepare the input for PCA
	double **cluster = new double*[SIZE];
	
	int level;
	if (i==-1 || j==-1) {
		level = (i==-1) ? levelBetaClusters[j] : levelBetaClusters[i];
	} else {
		level = (levelBetaClusters[i] > levelBetaClusters[j]) ? levelBetaClusters[i] : levelBetaClusters[j];
	}	
	
	Db *db = calcTree->getDb(level);
	
	// pointers to the parents of a cell
	stCell **parentsVector = new stCell*[level];
	for (int k=0;k<level;k++) {
		parentsVector[k] = new stCell();
	}//end for
	
	//prepare the fullId array
	int nPos = (int) ceil((double)DIM/8);
	unsigned char *fullId = new unsigned char[(level+1)*nPos];
	memset(fullId,0,(level+1)*nPos);
	
	//prepare the cell and the Dbts to receive 
	//key/data pairs from the dataset 
	stCell *cell = new stCell();
	stCellId *id = new stCellId();
	Dbt searchKey, searchData;
	searchKey.set_data(fullId);
	searchKey.set_ulen((level+1)*nPos);
	searchKey.set_flags(DB_DBT_USERMEM);
	searchData.set_data(cell);
	searchData.set_ulen(sizeof(stCell));
	searchData.set_flags(DB_DBT_USERMEM);
	
	// Get a cursor
	Dbc *cursorp;
	db->cursor(NULL,&cursorp,0);
	
	//prepare to walk through the level
	int belongsTo;
	double *normalizeSlope = calcTree->getNormalizeSlope(), *normalizeYInc = calcTree->getNormalizeYInc();
	double *maxCell = new double[DIM];
	double *minCell = new double[DIM];
	
	// iterate over the database, retrieving each record in turn
	int ret;
	while ((ret = cursorp->get(&searchKey, &searchData, DB_NEXT)) == 0) {
		//set id for cell
		id->setIndex(fullId+(level*nPos)); //copy from fullId to id
		cell->setId(id); //copy from id to cell->id
		//finds the parents of cell
		calcTree->findParents(fullId,parentsVector,level);
		// discovers the position of cell in the data space
		cellPosition(cell,parentsVector,minCell,maxCell,level);
		
		belongsTo=0;
		for (int betaCluster=0; (!belongsTo) && betaCluster < numBetaClusters; betaCluster++) {
			//test if this cluster is part of the PCA input
			if ( (betaCluster == i) || (betaCluster == j) 
				//|| (i != -1 && correlationClustersBelongings[i] != -1 && correlationClustersBelongings[betaCluster] == correlationClustersBelongings[i]) 
				//|| (j != -1 && correlationClustersBelongings[j] != -1 && correlationClustersBelongings[betaCluster] == correlationClustersBelongings[j]) 
				) {
				
				belongsTo=1;
				// verify if the cell belongs to the beta-cluster
				for (int dim=0; belongsTo && dim<DIM; dim++) {				
					if (! (  (((maxCell[dim]-minCell[dim])/2)+minCell[dim]) >= minBetaClusters[betaCluster][dim] && 
						     (((maxCell[dim]-minCell[dim])/2)+minCell[dim]) <= maxBetaClusters[betaCluster][dim] ) ) {
						belongsTo=0; // this cell does not belong to the beta-cluster
					}//end if
				}//end for
			}
		}
		
		if (belongsTo) { // this cell belongs to the PCA input
			for (int p=0; p < cell->getSumOfPoints(); p++) {
				cluster[*clusterSize] = new double[DIM];
				for (int dim=0; dim<DIM; dim++) {
					cluster[*clusterSize][dim] = (((((maxCell[dim]-minCell[dim])/2)+minCell[dim])*normalizeSlope[dim])+normalizeYInc[dim]);
				}
				(*clusterSize)++;
			}
		}
		
	}
	if (ret != DB_NOTFOUND) { //it should never enter here
		cout << "Error!" << endl;
	}
	
	//closes the cursor
	cursorp->close();
	
	//copy cluster to clusterMat (OpenCV format)
	cv::Mat clusterMat(*clusterSize,DIM,CV_64F);
	for (int p=0; p<*clusterSize; p++) {
		for (int d=0; d<DIM; d++) {
			clusterMat.at<double>(p,d) = cluster[p][d];
		}
		delete [] cluster[p]; //disposes the used memory
	}
	
	//disposes the used memory
	delete [] cluster;
	delete [] minCell;
	delete [] maxCell;
	delete cell;
	delete id;
	delete [] fullId;
	for (int k=0;k<level;k++) {
		delete parentsVector[k];
	}//end for
	delete [] parentsVector;
	
	return clusterMat;
	
}
