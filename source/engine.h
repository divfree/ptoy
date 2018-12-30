#pragma once

#include <cmath>
#include <vector>
#include <cassert>
#include <array>
#include <numeric>

#include "geom.h"

class Engine {
 public:
  void Step(Cont<Vect>& pp, Cont<Vect>& vv) {
    for (Size q = 0; q < pp.size(); ++q) {
      auto p = pp[q];
      vv[q] += (Vect(0.3, 0.) - p) * 0.01;
      pp[q] += vv[q];
    }
  }
};
