#include "Normalization.h"

#include "arboretum/stCommon.h"
#include <algorithm>
#include <iostream>

namespace Halite {
  template<typename D>
  Normalization<D>::Normalization(PointSource<D>& data, NormalizationMode mode) {
    size_t DIM=data.dimension();

    slope.resize(DIM, 1.0);
    yIntercept.resize(DIM, 0.0);
    
    std::vector<D> minD(DIM,0.0);
    std::vector<D> maxD(DIM,0.0);

    // normalizes the data
    minMax(data, minD, maxD);

    
    D normalizationFactor = 1.0;

    std::vector<D> resultPoint(DIM,0.0);
    for (size_t i = 0; i < DIM; i++) {
      slope[i] = (maxD[i] - minD[i]) * normalizationFactor; //a takes the range of each dimension
      yIntercept[i] = minD[i];
      if (slope[i] == 0) {
	slope[i] = 1;
      }//end if
    }//end for

    if (mode != Independent) {
      D slopeVal;
      if(mode==MaintainProportion) {
	slopeVal=*std::max_element(slope.begin(), slope.end());
      } else if(mode==Clip){
	slopeVal=*std::min_element(slope.begin(), slope.end());
      } else {
	std::cerr<<"Error! Unrecognized normalization mode\n";
      }
      slope.assign(slope.size(), slopeVal);
    }//end if
  
  }
  template<typename D>
  void Normalization<D>::minMax(PointSource<D>& data, std::vector<D>& min, std::vector<D>& max) {
    size_t DIM=data.dimension();
    for (size_t j = 0; j<DIM; j++){ // sets the values to the minimum/maximum possible here
      min[j] = std::numeric_limits<D>::max();
      max[j] = std::numeric_limits<D>::lowest();
    }// end for
    // looking for the minimum and maximum values
    for (data.restartIteration(); data.hasNext();) {
      const D* onePoint=data.readPoint();
      for (size_t j = 0; j<DIM; j++) {
	if (onePoint[j] < min[j]) {
	  min[j] = onePoint[j];
	}//end if
	if (onePoint[j] > max[j]) {
	  max[j] = onePoint[j];
	}//end if
      }//end for

    }//end for
  }

  template class Normalization<float>;
  template class Normalization<double>;
}
