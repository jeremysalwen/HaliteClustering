#include "Normalization.h"

#include "arboretum/stCommon.h"
#include <algorithm>
#include <iostream>

namespace Halite {
  Normalization::Normalization(PointSource& data, Mode mode) {
    size_t DIM=data.dimension();

    slope.resize(DIM, 1.0);
    yIntercept.resize(DIM, 0.0);
    
    std::vector<double> minD(DIM,0.0);
    std::vector<double> maxD(DIM,0.0);

    // normalizes the data
    minMax(data, minD, maxD);

    
    double normalizationFactor = 1.0;

    std::vector<double> resultPoint(DIM,0.0);
 
    for (size_t i = 0; i < DIM; i++) {
      slope[i] = (maxD[i] - minD[i]) * normalizationFactor; //a takes the range of each dimension
      yIntercept[i] = minD[i];
      if (slope[i] == 0) {
	yIntercept[i] = 1;
      }//end if
    }//end for

    if (mode != Independent) {
      double slopeVal;
      if(mode==MaintainProportion || mode == GeoReferenced) {
	slopeVal=*std::max_element(slope.begin(), slope.end());
      } else if(mode==Clip){
	slopeVal=*std::min_element(slope.begin(), slope.end());
      } else {
	std::cerr<<"Error! Unrecognized normalization mode\n";
      }
      slope.assign(slope.size(), slopeVal);
    }//end if
  
  }
   void Normalization::minMax(PointSource& data, std::vector<double>& min, std::vector<double>& max) {
    size_t DIM=data.dimension();
    for (size_t j = 0; j<DIM; j++){ // sets the values to the minimum/maximum possible here
      min[j] = MAXDOUBLE;
      max[j] = -MAXDOUBLE;
    }// end for
    // looking for the minimum and maximum values
    for (data.restartIteration(); data.hasNext();) {
      const double* onePoint=data.readPoint();
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
}
