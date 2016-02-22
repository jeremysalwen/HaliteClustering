#include "Classifier.h"

namespace Halite {
  template<typename D>
  void Classifier<D>::denormalize() {
    if(normalization) {
      for(BetaCluster<D>& b:betaClusters) {
	normalization->normalize(b.min.begin(), b.min.begin());
	normalization->normalize(b.max.begin(), b.max.begin());
      }
      normalization.reset();
    }
  }
  template class Classifier<float>;
  template class Classifier<double>;
}
