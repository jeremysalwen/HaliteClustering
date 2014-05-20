/**********************************************************************
 * GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
 *
 *                           Homepage: http://gbdi.icmc.usp.br/arboretum
 **********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * CiÍncias Matem·ticas e de ComputaÁ„o, University of S„o Paulo -
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
 *        de Dados e Imagens, Instituto de CiÍncias Matem·ticas e de
 *        ComputaÁ„o, University of S„o Paulo - Brazil (the Databases
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
 * This file demonstrates the usage of the stCorrelationClustering class.
 *
 * @version 1.0
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 */
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// MrCC
//------------------------------------------------------------------------------

#include <string.h>
#include <time.h>

#ifndef __GNUG__
#include "stCorrelationClustering.h"
#endif //__GNUG__

#ifdef __GNUG__
#include "stCorrelationClustering.cpp"
#endif //__GNUG__


// default values
#define TAM_LINE 500
#define NORMALIZE_FACTOR 0 // Independent

// global variables
clock_t startTime;
DBTYPE dbType = DB_HASH; //DB_HASH or DB_BTREE
char memory = 2; // 0: unlimited
                 // 1: limited
                 // 2: none

/**
 * Initiates the measurement of run time.
 */
void initClock() {
	startTime=clock();
}//end initClock

/**
 * Prints the elapsed time.
 */
void printElapsed() { 
	printf("Elapsed time: %0.3lf sec.\n\n",(double)(clock()-startTime)/CLOCKS_PER_SEC);
}//end printElapsed

/**
 * Initiates the clustering process.
 *
 * @param argc Number of arguments received.
 * @param argv Array with the arguments received.
 * @return success (0), error (1).
 */
int main(int argc, char **argv) {
	initClock(); // initiates the meassurement of run time
	
	// first validations
	if (argc != 5) {
		printf("Usage: MrCC <pThreshold> <H> <hardClustering> <initialLevel>");
		return 1; //error
	}//end if
	
	if (atoi(argv[2]) < 2) {
		printf("MrCC needs at least two resolution levels (H >= 2) to perform the clustering process.");
		return 1; //error
	}//end if
	
	// opens/creates the used files
	FILE *database, *result;
    database=fopen(INPUT, "r");
	result=fopen(OUTPUT, "w");
	if (!(database&&result)) {
		printf("MrCC could not open the database file or create the result file.");
		return 1; //error
	}//end if
	
	double **objectsArray = 0;
	if (memory == 0) { //unlimited memory
		// reads objects from the source database
		objectsArray = new double*[SIZE];
		for (int i=0; i<SIZE; i++) {
			objectsArray[i] = new double[DIM];
			readPoint(database, objectsArray[i]);
		}//end for
	}

	// creates an object of the class stCorrelationClustering
    stCorrelationClustering *sCluster = new stCorrelationClustering(objectsArray, database, NORMALIZE_FACTOR, 
										   (2*DIM), -1, atof(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), dbType, memory);		
	
	printf("The tree was built.\n");
	printElapsed(); // prints the elapsed time
	
	printf("Time spent in the normalization.\n");
	printf("Time: %0.3lf sec.\n\n",(double)(sCluster->getTimeNormalization())/CLOCKS_PER_SEC);
	
    // looks for correlation clusters
    sCluster->findCorrelationClusters();
	
	printf("The clusters were found.\n");
	printElapsed(); // prints the elapsed time
	
	// mounts the result file
	int numBetaClusters = sCluster->getNumBetaClusters(), numCorrelationClusters = sCluster->getNumCorrelationClusters(), betaCluster,
	*correlationClustersBelongings = sCluster->getCorrelationClustersBelongings();
	char **dimCorrelationClusters = sCluster->getDimCorrelationClusters(), line[TAM_LINE], belongsTo, str[20];
	double **minBetaClusters = sCluster->getMinBetaClusters(), **maxBetaClusters = sCluster->getMaxBetaClusters(),
	*normalizeSlope = sCluster->getCalcTree()->getNormalizeSlope(),
	*normalizeYInc = sCluster->getCalcTree()->getNormalizeYInc();
	
	// axes relevant to the found clusters
	for (int i=0; i<numCorrelationClusters; i++) {
		strcpy(line,""); // empty line
		strcat(line,"ClusterResult");
		sprintf(str,"%d",i+1);		
		strcat(line,str);		
		for (int j=0; j<DIM; j++) {
			(dimCorrelationClusters[i][j]) ? strcat(line," 1") : strcat(line," 0");
		}//end for
		strcat(line,"\n");
		fputs(line,result); // writes the relevant axes to the current cluster in the result file
	}//end for
	
	// labels each point after the clusters found
	fputs("LABELING\n",result);
	double *onePoint = new double[DIM];
	if (memory != 0) { // limited or none memory
		fseek(database,0,SEEK_SET); //go to the file begin
	}
	if (atoi(argv[3])) { //hard clustering
		for (int point=0; point < SIZE; point++) {
			(memory == 0) ? copyPoint(objectsArray[point], onePoint) : readPoint(database, onePoint);
			strcpy(line,""); // empty line
			belongsTo=0;
			for (betaCluster=0; (!belongsTo) && betaCluster < numBetaClusters; betaCluster++) {
				belongsTo=1;
				// undoes the normalization and verify if the current point belongs to the current beta-cluster
				for (int dim=0; belongsTo && dim<DIM; dim++) {				
					if (! (onePoint[dim] >= ((minBetaClusters[betaCluster][dim]*normalizeSlope[dim])+normalizeYInc[dim]) && 
						   onePoint[dim] <= ((maxBetaClusters[betaCluster][dim]*normalizeSlope[dim])+normalizeYInc[dim])) ) {
						belongsTo=0; // this point does not belong to the current beta-cluster
					}//end if
				}//end for
			}//end for
			if (belongsTo) { // this point belongs to one cluster
				sprintf(line,"%d %d\n",point+1, correlationClustersBelongings[betaCluster-1]+1);			
			} else {
				sprintf(line,"%d %d\n",point+1,0);
			}//end if
			fputs(line,result);
		}//end for
	} else { //soft clustering
		int outlier;
		for (int point=0; point < SIZE; point++) {
			(memory == 0) ? copyPoint(objectsArray[point], onePoint) : readPoint(database, onePoint);
			outlier = 1;
			for (betaCluster=0; betaCluster < numBetaClusters; betaCluster++) {
				belongsTo=1;
				// undoes the normalization and verify if the current point belongs to the current beta-cluster
				for (int dim=0; belongsTo && dim<DIM; dim++) {				
					if (! (onePoint[dim] >= ((minBetaClusters[betaCluster][dim]*normalizeSlope[dim])+normalizeYInc[dim]) && 
						   onePoint[dim] <= ((maxBetaClusters[betaCluster][dim]*normalizeSlope[dim])+normalizeYInc[dim])) ) {
						belongsTo=0; // this point does not belong to the current beta-cluster
					}//end if
				} // end for
				if (belongsTo) {
					outlier = 0;
					sprintf(str, "%d %d\n", point+1, correlationClustersBelongings[betaCluster]+1);
					fputs(str, result);
				}//end if
			}//end for
			if (outlier) {
				sprintf(str, "%d %d\n", point+1, 0);
				fputs(str, result);
			}
		}//end for
	}
	fclose(result); // the result file will not be used anymore
	fclose(database); // the database file will not be used anymore
	
	printf("The result file was built.\n");
	printElapsed(); // prints the elapsed time

	// disposes the used structures
	delete [] onePoint;
	if (memory == 0) { //unlimited memory
		// disposes objectsArray
		for (int i=0;i<SIZE;i++) {
			delete [] objectsArray[i];
		}//end for
		delete [] objectsArray;		
	}
	delete sCluster;
	
	printf("The used structures were deleted. Done!\n");
	printElapsed(); // prints the elapsed time
	//getc(stdin);
	return 0; // success
}//end main
