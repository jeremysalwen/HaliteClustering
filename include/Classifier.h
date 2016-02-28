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
      if(betaClusters.empty()) {
	return out;
      }      
      size_t DIM=betaClusters[0].min.size();

      std::vector<D> normalized(normalization?DIM:0);
      if(normalization) {
	normalization->normalize(point, normalized.data());
      }
      
      for(size_t i=0; i<betaClusters.size(); i++) {
	BetaCluster<D>& betaCluster=betaClusters[i];

	bool contains;
	if(normalization) {
	  contains=betaCluster.contains(normalized.begin());
	} else {
	  contains=betaCluster.contains(point);
	}
	
	if(contains) {
	  *out++ = betaCluster.correlationCluster+1;
	  if(this->hardClustering) {
	    return out;
	  }
	}	
      }
      return out;
    }

    void denormalize();
    
    bool hardClustering;

    std::shared_ptr<Normalization<D>> normalization;

    std::vector<BetaCluster<D> > betaClusters;
    };
}

#endif
