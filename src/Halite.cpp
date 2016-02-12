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
 * This file demonstrates the usage of the haliteClustering class.
 *
 * @version 1.0
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 */
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// Halite
//------------------------------------------------------------------------------

#include <string.h>
#include <time.h>

#include "haliteClustering.h"

#include "PointSource.h"

using namespace Halite;

#define INPUT "databases/ds2.csv" //input data path
#define OUTPUT "results/resultds12.dat" //output path

// default values
#define NORMALIZE_FACTOR 0 // Independent

// global variables
clock_t startTime;
DBTYPE dbType = DB_HASH; //DB_HASH or DB_BTREE

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
	if (argc != 6) {
		printf("Usage: Halite <pThreshold> <H> <hardClustering> <initialLevel> <memoryMode>\n");
		return 1; //error
	}//end if
	
	if (atoi(argv[2]) < 2) {
		printf("Halite needs at least two resolution levels (H >= 2) to perform the clustering process.\n");
		return 1; //error
	}//end if
	
	char memory=atoi(argv[5]);
	if (memory<0 || memory >3) {
		printf("Possible memory modes are 0: Everything in memory, 1: Only Data in Memory 2: Only Tree in Memory 3: Everything on disk\n");
		return 1;
	}//end if

	// opens/creates the used files
	FILE  *result;
	result=fopen(OUTPUT, "w");
	if (!result) {
		printf("Halite could not create the result file.\n");
		return 1; //error
	}//end if

	PointSource*  db;
	try {
		db=new TextFilePointSource(INPUT);
	} catch(std::exception& e) {
		printf("'Halite could not open database file.\n");
		return 1;
	}
	size_t DIM=db->dimension();
	PointSource* memdb=NULL;
	PointSource* datasource=db;
	std::vector<double*> objectsArray;
	if (!(memory & 1)) { //unlimited memory
		for(db->restartIteration(); db->hasNext();) {
			double* tmp = new double[DIM];
			const double* t=db->readPoint();
			for(unsigned int j=0; j<DIM; j++) {
				tmp[j]=t[j];
			}
			objectsArray.push_back(tmp);
		}//end for
		memdb=new ArrayOfPointersPointSource(&objectsArray[0],db->dimension(),objectsArray.size());
		datasource=memdb;
	}

	// creates an object of the class haliteClustering
    haliteClustering *sCluster = new haliteClustering(*datasource, NORMALIZE_FACTOR, (2*DIM), -1, atof(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), dbType, (memory & 2));		
	
	printf("The tree was built.\n");
	printElapsed(); // prints the elapsed time
	
	printf("Time spent in the normalization.\n");
	printf("Time: %0.3lf sec.\n\n",(double)(sCluster->getTimeNormalization())/CLOCKS_PER_SEC);
	
    // looks for correlation clusters
    sCluster->findCorrelationClusters();
	
	printf("The clusters were found.\n");
	printElapsed(); // prints the elapsed time
	
	// mounts the result file

	unsigned int numCorrelationClusters=sCluster->getNumCorrelationClusters();

        std::vector<char*> dimCorrelationClusters = sCluster->getDimCorrelationClusters();
	
	// axes relevant to the found clusters
	for (unsigned int i=0; i<numCorrelationClusters; i++) {
		fputs("ClusterResult",result);
		fprintf(result,"%d",i+1);	
		for (unsigned int j=0; j<DIM; j++) {
			(dimCorrelationClusters[i][j]) ? fputs(" 1",result) : fputs(" 0",result);
		}//end for
		fputs("\n",result); // writes the relevant axes to the current cluster in the result file
	}//end for
	
	// labels each point after the clusters found
	fputs("LABELING\n",result);
	const double *onePoint;

	int point=0;
	for (datasource->restartIteration(); datasource->hasNext(); point++) {
		onePoint=datasource->readPoint();
		std::vector<int> clusters;
		sCluster->assignToClusters(onePoint,std::back_inserter(clusters));
		for(unsigned int i=0; i<clusters.size(); i++) {
			fprintf(result, "%d %d\n", point+1, clusters[i]);
		}
		if(clusters.empty()) {
			fprintf(result, "%d %d\n", point+1, 0);
		}
	}//end for

	fclose(result); // the result file will not be used anymore
	delete db;

	printf("The result file was built.\n");
	printElapsed(); // prints the elapsed time

	// disposes the used structures
	if (!(memory & 1)) { //unlimited memory
		// disposes objectsArray
		for (unsigned int i=0;i<objectsArray.size();i++) {
			delete [] objectsArray[i];
		}//end for	
		delete memdb;	
	}
	delete sCluster;
	
	printf("The used structures were deleted. Done!\n");
	printElapsed(); // prints the elapsed time
	//getc(stdin);
	return 0; // success
}//end main
