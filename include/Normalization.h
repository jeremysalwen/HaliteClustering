//-*-c++-*-
#ifndef NORMALIZATION_H
#define NORMALIZATION_H

#include <vector>
#include "PointSource.h"


namespace Halite {
  
  class Normalization {
  public:
    /* when we have the proportional normalization, every dimension is gonna receive
       the biggest range.
       In Clip mode, they all receive the smallest range.
       In Independent mode, they are all independently sized
    */
    enum Mode {
      Independent = 0, MaintainProportion, Clip, GeoReferenced
    };
    Normalization(size_t DIM) : slope(DIM,1), yIntercept(DIM,0){}
    Normalization(PointSource& data, Mode mode);

    template<typename it, typename it2>
    void normalize(it in, it2 out) {
      for(size_t i=0; i<slope.size(); i++) {
	out[i]=(in[i]-yIntercept[i])/slope[i];
      }
    }
    
    template<typename it, typename it2>
    void denormalize(const it in, it2 out) {
      for(size_t i=0; i<slope.size(); i++) {
	out[i]=in[i]*slope[i] + yIntercept[i];
      }
    }
    
    /**
     * Finds the minimum and maximum data values in each dimension.
     * Method from the stFractalDimension class (adapted).
     *
     * @param PointSource Source of the database objects.
     * @param min Vector to receive the minimum data value in each dimension.
     * @param max Vector to receive the maximum data value in each dimension.
     *
     */

    static void minMax(PointSource& data, std::vector<double>& min, std::vector<double>& max);
    std::vector<double> slope;
    std::vector<double> yIntercept;
  };
}
#endif
