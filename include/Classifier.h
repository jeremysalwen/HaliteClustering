// -*-c++-*-
#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include <vector>
#include <cstddef>
#include <memory>
#include "BetaCluster.h"
#include "Normalization.h"

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
      Iterator assignToClusters(const D* point, Iterator out) {
      for(size_t i=0; i<betaClusters.size(); i++) {
	BetaCluster<D>& betaCluster=betaClusters[i];
	bool belongsTo=true;
	std::vector<D> denormMin(betaCluster.min.size());
	std::vector<D> denormMax(betaCluster.max.size());
	normalization->denormalize(betaCluster.min.begin(), denormMin.begin());
	normalization->denormalize(betaCluster.max.begin(), denormMax.begin());
	
	// undoes the normalization and verify if the current point belongs to the current beta-cluster
	for (size_t dim=0; belongsTo && dim<betaCluster.min.size(); dim++) {			       
	  if (! (point[dim] >= denormMin[dim] && 
		 point[dim] <= denormMax[dim]) ) {
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

    std::shared_ptr<Normalization<D>> normalization;

    std::vector<BetaCluster<D> > betaClusters;
  };
}

#endif
