// -*-c++-*-
#ifndef CORRELATION_CLUSTER_H
#define CORRELATION_CLUSTER_H

#include <vector>

namespace Halite {
  class CorrelationCluster {
  public:
  CorrelationCluster(int DIM) :relevantDimension(DIM, false) { 
    }
    //vector<bool>, except not specialized
    std::vector<unsigned char> relevantDimension;
  };
}

#endif
