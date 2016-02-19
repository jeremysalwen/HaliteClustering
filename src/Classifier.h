// -*-c++-*-
#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include <vector>
#include <cstddef>
#include "BetaCluster.h"
namespace Halite {

  template <typename D>
  class Classifier {
  public:
    /**
     * Finds clustering of a point, and writes all clusters to the out iterator
     * if you are using hard clustering, just pass a pointer to the integer you want it to write to (it will not write anything for outliers)
     * if you are using soft clustering, pass a std::back_inserter to your favorite container to store the results.
     */
    template<typename Iterator>
      Iterator assignToClusters(const double* point, Iterator out) {
      for(size_t i=0; i<betaClusters.size(); i++) {
	BetaCluster<double>& betaCluster=betaClusters[i];
	bool belongsTo=true;
	// undoes the normalization and verify if the current point belongs to the current beta-cluster
	for (size_t dim=0; belongsTo && dim<betaCluster.min.size(); dim++) {			       
	  if (! (point[dim] >= ((betaCluster.min[dim]*normalizeSlope[dim])+normalizeYInc[dim]) && 
		 point[dim] <= ((betaCluster.max[dim]*normalizeSlope[dim])+normalizeYInc[dim])) ) {
	    belongsTo = false; // this point does not belong to the current beta-cluster
	  }
	}
	if(belongsTo) {
	  *out++ = betaCluster.correlationCluster+1;
	  if(this->hardClustering) {
	    return out;
	  }
	}	
      }
      return out;
    }

    
    bool hardClustering;

    const double* normalizeSlope;
    const double* normalizeYInc;

    std::vector<BetaCluster<D> > betaClusters;
  };
}

#endif
