// -*-c++-*-
#ifndef BETA_CLUSTER_H
#define BETA_CLUSTER_H

#include <vector>
#include <cstddef>
namespace Halite {
  template <typename D>
    class BetaCluster {
  public:
  BetaCluster(int levl, int DIM): level(levl), cost(-1), correlationCluster(-1), relevantDimension(DIM, false), min(DIM), max(DIM) {
    }

    bool contains(std::vector<D>& point) const {
      for(size_t i=0; i<point.size(); i++){
	if(point[i]<min[i] || point[i]>max[i]) {
	  return false;
	}
      }
      return true;
    }
    int level;
    int cost;
    int correlationCluster;
    //vector<bool>, except not specialized
    std::vector<unsigned char> relevantDimension;
    std::vector<D> min;
    std::vector<D> max;
  };
}
#endif
