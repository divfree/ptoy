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
    auto ppt = pp;
    auto vvt = vv;
    {
      auto ff = Rhs(pp);
      for (Size q = 0; q < pp.size(); ++q) {
        Scal dt = kTimeStep;
        vvt[q] += ff[q] * dt;
        ppt[q] += vvt[q] * dt;
      }
    }

    {
      auto ff = Rhs(ppt);
      for (Size q = 0; q < pp.size(); ++q) {
        Scal dt = kTimeStep;
        vv[q] += ff[q] * dt;
        pp[q] += vv[q] * dt;
      }
    }
  }
  Vect F12(Vect p1, Vect p2, Scal R, Scal sigma) {
    const Vect dp = p1 - p2;
    const Scal r2 = dp.dot(dp);
    const Scal r2inv = 1. / r2;
    const Scal d2 = r2inv * (R * R);
    const Scal d6 = d2 * d2 * d2;
    const Scal d12 = d6 * d6;
    return dp * std::max<Scal>(0., sigma * (d12 - d6) * r2inv);
  }
  // pp: target positions
  // n: target number
  // ppo: other positions
  // no: other number
  void CalcForceSerial(Vect* ff, const Vect* pp, Size n, 
                       const Vect* ppo, Size no) {
    for (Size i = 0; i < n; ++i) {
      for (Size j = 0; j < no; ++j) {
        if (&pp[i] != &ppo[j]) {
          ff[i] += F12(pp[i], ppo[j], kRadius, kSigma);
        }
      }
    }
  }
  Cont<Vect> Rhs(const Cont<Vect>& pp) {
    Cont<Vect> ff(pp.size(), Vect(0));
    CalcForceSerial(ff.data(), pp.data(), pp.size(), pp.data(), pp.size());

    for (Size q = 0; q < pp.size(); ++q) {
      Vect r = -pp[q];
      if (r.norm() > kRadius) {
        ff[q] -= r * (-kPointForceAttractive / std::pow(r.norm(), 3));
      }
    }

    return ff;
  }
};
