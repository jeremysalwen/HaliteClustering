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
* This file defines the class stCorrelationClustering.
*
* @version 1.0
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Caetano Traina Jr (caetano@icmc.usp.br)
*
*/
// Copyright (c) 2002-2009 GBDI-ICMC-USP

#ifndef __STCORRELATIONCLUSTERING_H
#define __STCORRELATIONCLUSTERING_H

#include <cv.h> //OpenCV
#include "arboretum/stCountingTree.h"

#ifndef __GNUG__
	#include "Utile.h"	
#endif //__GNUG__

#ifdef __GNUG__
	#include "Utile.cpp"
#endif //__GNUG__

#include <math.h>
#include <stdio.h>
#include <time.h>

//----------------------------------------------------------------------------
// class stCorrelationClustering
//----------------------------------------------------------------------------
/**
* This class is used to find clusters in subspaces of the original data space.
*
* @version 1.0
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Caetano Traina Jr (caetano@icmc.usp.br)
*/
//---------------------------------------------------------------------------
class stCorrelationClustering {

   public:

      /**
      * Creates the needed structures to find correlation clusters.
      *
      * @param objectsArray Array with the database objects.
      * @param normalizeFactor Determines how data will be normalized.
      * @param centralConvolutionValue Determines the central weight in the convolution matrix.
      * @param neighbourhoodConvolutionValue Determines the face neighbours weight in the convolution matrix.
      * @param pThreshold Threshold used to spot a beta-cluster, based on the binomial probability.
      * @param H The number of grid levels to build the counting tree.
	  * @param hardClustering Choose between hard (1) and soft (0) clustering.
      *
      */
      stCorrelationClustering (double ** objectsArray, FILE *database,
							   int normalizeFactor,
                               int centralConvolutionValue, int neighbourhoodConvolutionValue,
                               double pThreshold, int H, int hardClustering, int initialLevel, DBTYPE dbType, char memory);

      /**
      * Disposes the allocated memory.
      */
      ~stCorrelationClustering();

      /**
      * Finds clusters in subspaces.
      */
      void findCorrelationClusters();

      /**
      * Gets the number of beta-clusters found.
      */
      int getNumBetaClusters() {
        return numBetaClusters;
      }//end getNumBetaClusters

      /**
      * Gets the number of correlation clusters (merged).
      */
      int getNumCorrelationClusters() {
        return numCorrelationClusters;
      }//end getNumCorrelationClusters

      /**
      * Gets the minimum bounds of the beta-clusters found.
      */
      double **getMinBetaClusters() {
        return minBetaClusters;
      }//end getMinBetaClusters

      /**
      * Gets the maximum bounds of the beta-clusters found.
      */
      double **getMaxBetaClusters() {
        return maxBetaClusters;
      }//end getMaxBetaClusters

      /**
      * Gets the relevant dimensions to the correlation clusters.
      */
      char **getDimCorrelationClusters() {
        return dimCorrelationClusters;
      }//end getDimCorrelationClusters

      /**
      * Gets the correlation clusters belongings.
      */
      int *getCorrelationClustersBelongings() {
        return correlationClustersBelongings;
      }//end getCorrelationClustersBelongings

      /**
      * Gets the used counting tree.
      */
      stCountingTree *getCalcTree() {
        return calcTree;
      }//end getCalcTree
	
	  /**
	  * Gets the time spent in the normalization.
	  */
	  clock_t getTimeNormalization() {
		return timeNormalization;
	  }//end getTimeNormalization

   private:
	
	  /**
	  * Time spent in the normalization.
	  */
	  clock_t timeNormalization;

      /**
      * Counting-tree pointer.
      */
      stCountingTree *calcTree;


      /**
      * Number of grid levels in the counting tree.
      */
      int H;

      /**
      * Choose between hard and soft clustering.
      */
      int hardClustering;
	
	  /**
	  * Defines the initial tree level to look for clusters.
	  */
	  int initialLevel;

      /**
      * Pointers to the central cell of a found beta-cluster and to a neighbour cell.
      */
      stCell *betaClusterCenter;
      stCell *neighbour;

      /**
      * Vectors to describe the positions of the beta-cluster center in the data space.
      */
      double *minBetaClusterCenter;
      double *maxBetaClusterCenter;
      
	  /**
      * Position of a neighbour of the beta-cluster center in the data space, regarding a dimension e_j.
      */      
	  double minNeighbour;
      double maxNeighbour;

      /**
      * Vectors with pointers to the parents of a beta-cluster center and of a neighbour.
      */
      stCell **betaClusterCenterParents;
      stCell **neighbourParents;

      /**
      * Vectors used when discovering relevant attributes.
      */
      double *attributesRelevance;      

      /**
      * Vector used to indicate which direct neighbours of a central cell also belong to the beta-cluster.
      */
      char *neighbourhood;

      /**
      * Defines the convolution matrix (center and direct neighbours).
      */
      int centralConvolutionValue;
      int neighbourhoodConvolutionValue;

      /**
      * Defines the threshold used to spot a beta-cluster, based on the binomial probability.
      */
      double pThreshold;

      /**
      * Number of beta-clusters found.
      */
      int numBetaClusters;

      /**
      * Number of correlation clusters (merged).
      */
      int numCorrelationClusters;

      /**
      * Pointers to the vectors that describe the beta-clusters and the correlation clusters found.
      */
      double **minBetaClusters; // array L in the paper
      double **maxBetaClusters; // array U in the paper
      char **dimBetaClusters; // array V in the paper
	  char **dimCorrelationClusters;
	  int *correlationClustersBelongings;

	  /**
      * Merges beta-clusters that share some database space.
      */
      void mergeBetaClusters();
	
	  int shouldMerge(int i, int j);
	
	  int cost(int i, int j);
	
	  int indCost(double n);
	
	  cv::Mat inputPCA(int i, int j, int *clusterSize);
	
	  int *costBetaClusters;
	  int *levelBetaClusters;

      /**
      * Calculates the cThreshold based on the Minimun Description Length (MDL) method.
      *
      * @param attributesRelevance Vector with the calculed relevances of each attribute.
      *
      */
      double calcCThreshold(double *attributesRelevance);

      /**
      * Finds the best cut point position based on the MDL method.
      *
      * @param sortedRelevance Vector with the sorted relevances of each attribute.
      *
      */
      int minimumDescriptionLength(double *sortedRelevance);

      /**
      * Walk through the counting tree applying the convolution matrix
      * to each cell in a defined level.
      *
      * @param level The counting tree level to be analyzed.
      *
      */
      int walkThroughConvolution(int level);

      /**
      * Applies the convolution matrix to a grid cell.
      *
      * @param cell The grid cell to apply the convolution matrix.
      * @param cellParents The cell parents.
      * @param level The cell level in the counting tree.
      *
      */
      int applyConvolution(stCell *cell, stCell **cellParents,
                                          int level);

      /**
      * Finds the position of a cell in the data space.
      *
      * @param cell The analyzed cell.
      * @param cellParents The cell parents.
      * @param min Vector to receive the minimum position limits of cell in each dimension.
      * @param max Vector to receive the maximum position limits of cell in each dimension.
      * @param level The level of cell in the counting tree.
      *
      */
      void cellPosition(stCell *cell, stCell **cellParents,
                                        double *min, double *max, int level);

      /**
      * Finds the position of a cell in the data space, regarding a dimension e_j.
      *
      * @param cell The analyzed cell.
      * @param cellParents The cell parents.
      * @param min The minimum position limit of cell in dimension e_j.
      * @param max The maximum position limit of cell in dimension e_j.
      * @param level The level of cell in the counting tree.
	  * @param j The dimension to be considered.
      *
      */
      void cellPositionDimensionE_j(stCell *cell, stCell **cellParents,
                                        double *min, double *max, int level, int j);

      /**
      * Finds the external face neighbour of cell in a determined dimension.
      *
      * @param dimIndex The dimension to look for the neighbour.
      * @param cell The analyzed cell.
      * @param neighbour Pointer to to the memory space where the found neighbour is returned.	
      * @param cellParents Vector with the cell parents.
      * @param neighbourParents Vector to receive the found neighbour parents.
      * @param level The level of cell in the counting tree.
      *
	  * ps.: externalNeighbour assumes that cellParents and neighbourParents have the very
	  *      same containt when the function is called.
      */
      int externalNeighbour(int dimIndex, stCell *cell, stCell **neighbour,
                                                                  stCell **cellParents,
                                                                  stCell **neighbourParents, int level);

      /**
      * Finds the internal face neighbour of cell in a determined dimension.
      *
      * @param dimIndex The dimension to look for the neighbour.
      * @param cell The analyzed cell.
      * @param neighbour Pointer to to the memory space where the found neighbour is returned.
      * @param cellParents Vector with the cell parents.	   
      * @param level The level of cell in the counting tree.
      *
      */
      int internalNeighbour(int dimIndex, stCell *cell, stCell **neighbour, stCell **cellParents, int level);
      
	  /**
      * compare function used by the qsort.
      *
      * @param a first double element.
      * @param b second double element.
      *
      */
      static int compare (const void *a, const void *b) {	     
		  if(*(double*)a > *(double*)b) {
			 return 1;
		  } else {	
			  if(*(double*)a < *(double*)b) {
				  return -1;
			  }//end if
		  }//end if
		  return 0;
	  }//end compare

      /**
      * Normalize data points and insert them in the counting tree.
      * Method from the stFractalDimension class (adapted).
      *
      * @param objectsArray Array with the database objects.
      * @param normalizeFactor Determines how data will be normalized.
      *
      */
      void fastDistExponent(double **objectsArray, FILE *database, int normalizeFactor, char memory);

      /**
      * Finds the minimum and maximum data values in each dimension.
      * Method from the stFractalDimension class (adapted).
      *
      * @param objectsArray Array with the database objects.
      * @param min Vector to receive the minimum data value in each dimension.
      * @param max Vector to receive the maximum data value in each dimension.
      *
      */
      void minMax(double **objectsArray, FILE *database, double *min, double *max, char memory);

};//end stCorrelationClustering
//----------------------------------------------------------------------------

#endif //__STCORRELATIONCLUSTERING_H
