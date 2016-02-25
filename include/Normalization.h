//-*-c++-*-
#ifndef NORMALIZATION_H
#define NORMALIZATION_H

#include <vector>
#include "PointSource.h"


namespace Halite {

  /* when we have the proportional normalization, every dimension is gonna receive
       the biggest range.
       In Clip mode, they all receive the smallest range.
       In Independent mode, they are all independently sized
    */
    enum NormalizationMode {
      Independent = 0, MaintainProportion, Clip
    };

  
  template <typename D>
  class Normalization {
  public:

    Normalization(size_t DIM) : slope(DIM,1), yIntercept(DIM,0){}
    Normalization(PointSource<D>& data, NormalizationMode mode);

    void normalize(const D* in, D* out) {
      for(size_t i=0; i<slope.size(); i++) {
	out[i]=(in[i]-yIntercept[i])/slope[i];
      }
    }
    void normalize(typename std::vector<D>::const_iterator in, typename std::vector<D>::iterator out) {
      for(size_t i=0; i<slope.size(); i++) {
	out[i]=(in[i]-yIntercept[i])/slope[i];
      }
    }
    
    void denormalize(const D* in, D* out) {
      for(size_t i=0; i<slope.size(); i++) {
	out[i]=in[i]*slope[i] + yIntercept[i];
      }
    }

    void denormalize(typename std::vector<D>::const_iterator in,typename std::vector<D>::iterator out) {
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

    static void minMax(PointSource<D>& data, std::vector<D>& min, std::vector<D>& max);
    std::vector<D> slope;
    std::vector<D> yIntercept;
  };
}
#endif
