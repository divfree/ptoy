#pragma once

#include <cmath>
#include <vector>
#include <cassert>
#include <array>
#include <numeric>
#include <x86intrin.h>

#include "geom.h"
#include "cell.h"

class Engine {
 public:
  void SetForce(Vect p, bool s) {
    fr_p_ = p;
    fr_s_ = s;
  }
  void SetForce(bool s) {
    fr_s_ = s;
  }
  void SetForce(Vect p) {
    fr_p_ = p;
  }
  void Step(Cont<Vect>& pp, Cont<Vect>& vv, CellList& cl, Scal& t) {
    auto ppt = pp;
    auto vvt = vv;
    Scal dt = kTimeStep;

    {
      auto ff = Rhs(pp, vv, cl);
      Scal dth = dt * 0.5;
      for (Size q = 0; q < pp.size(); ++q) {
        vvt[q] += ff[q] * dth / kMass;
        ppt[q] += vvt[q] * dth;
      }
    }

    {
      auto ff = Rhs(ppt, vvt, cl);
      for (Size q = 0; q < pp.size(); ++q) {
        vv[q] += ff[q] * dt / kMass;
        pp[q] += vv[q] * dt;
      }
      t += dt;
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
                       const Vect* ppo, Size no) { for (Size i = 0; i < n; ++i) {
      for (Size j = 0; j < no; ++j) {
        if (pp[i] != ppo[j]) {
          ff[i] += F12(pp[i], ppo[j], kRadius, kSigma);
        }
      }
    }
  }
  void CalcForceAvx(Vect* ff, const Vect* pp, Size n, 
                    const Vect* ppo, Size no) {
    return;
    using M = const __m256;
    // sigma = kSigma;
    M sigma = _mm256_broadcast_ss(&kSigma);
    // R2 = (2. * kRadius) ^ 2;
    const float tmp = std::pow(2. * kRadius, 2);
    M R2 = _mm256_broadcast_ss(&tmp);
    // threshold = kRadius ^ 2 * 1e-3
    const float tmp_th = std::pow(kRadius, 2) * 1e-3;
    M threshold = _mm256_broadcast_ss(&tmp_th);
    const float tmp_zero = 0.;
    M zero = _mm256_broadcast_ss(&tmp_zero);

    for (Size j = 0; j < no; ++j) {
      M qx = _mm256_broadcast_ss((float*)&ppo[j][0]);
      M qy = _mm256_broadcast_ss((float*)&ppo[j][1]);
      // qxy = (q.x, q.y)
      M qxy = _mm256_blend_ps(qx, qy, 0xAA);
      for (Size i = 0; i < n; i += 8) {
        // pxy =(p.x, p.y)
        M pxy_l = _mm256_load_ps((float*)&pp[i]);
        M pxy_h = _mm256_load_ps((float*)&pp[i + 4]);
        // rxy = pxy - qxy
        M rxy_l = _mm256_sub_ps(pxy_l, qxy);
        M rxy_h = _mm256_sub_ps(pxy_h, qxy);
        // rxy2 = rxy * rxy 
        M rxy2_l = _mm256_mul_ps(rxy_l, rxy_l);
        M rxy2_h = _mm256_mul_ps(rxy_h, rxy_h);
        // r2 = (rx * rx + ry * ry)
        // r2 = ([7] [6] [3] [2] [5] [4] [1] [0]) 
        __m256 r2 = _mm256_hadd_ps(rxy2_l, rxy2_h);
        if (true) {
          // r2 = max(r2, threshold)
          r2 = _mm256_max_ps(r2, threshold);
        }
        // c2 = 1. / r2
        M c2 = _mm256_rcp_ps(r2);
        // d2 = c2 * R2 = R2 / r2 
        M d2 = _mm256_mul_ps(c2, R2);
        // d6 = d2 * d2 * d2
        M d6 = _mm256_mul_ps(d2, _mm256_mul_ps(d2, d2));
        // d12 = d6 * d6
        M d12 = _mm256_mul_ps(d6, d6);
        // k = (d12 - d6) * sigma * c2
        __m256 k = _mm256_mul_ps(
            sigma, _mm256_mul_ps(c2, _mm256_sub_ps(d12, d6)));

        k = _mm256_max_ps(k, zero);

        // lo = k([3] [3] [2] [2] [1] [1] [0] [0])
        M kxy_l = _mm256_unpacklo_ps(k, k);
        // hi = k([7] [7] [6] [6] [5] [5] [4] [4]) 
        M kxy_h = _mm256_unpackhi_ps(k, k);

        // load force to fxy
        __m256 fxy_l = _mm256_load_ps((float*)&ff[i]);
        __m256 fxy_h = _mm256_load_ps((float*)&ff[i + 4]);
        // fxy += rxy * kxy
        fxy_l =  _mm256_add_ps(fxy_l, _mm256_mul_ps(kxy_l, rxy_l));
        fxy_h =  _mm256_add_ps(fxy_h, _mm256_mul_ps(kxy_h, rxy_h));
        // store force
        _mm256_store_ps((float*)&ff[i], fxy_l);
        _mm256_store_ps((float*)&ff[i + 4], fxy_h);
      }
    }
  }
  Cont<Vect> Rhs(const Cont<Vect>& pp, const Cont<Vect>& vv, CellList& cl) {
    Cont<Vect> ff(pp.size(), Vect(0));
    const Vect d(1.);
    cl.Reinit(pp, MIdx(2. * d[0] / kBlockSize + 0.5, 
                       2. * d[1] / kBlockSize + 0.5), -d, d);
    const MIdx dims = cl.GetDims();
    auto& kkc = cl.GetOffset();
    auto& qq = cl.GetPartIdx();

    #pragma omp parallel for schedule(dynamic)
    for (Size j = 0; j < dims[1]; ++j) {
      for (Size i = 0; i < dims[0]; ++i) {
        Size c = j * dims[0] + i;

        for (int dj = -1; dj <= 1; ++dj) {
          for (int di = -1; di <= 1; ++di) {
            int ii = i + di;
            int jj = j + dj;
            if (ii >= 0 && ii < dims[0] && jj >= 0 && jj < dims[1]) {
              Size co = jj * dims[0] + ii;
              for (Size k = kkc[c]; k < kkc[c + 1]; ++k) {
                for (Size ko = kkc[co]; ko < kkc[co + 1]; ++ko) {
                  if (k != ko) {
                    ff[qq[k]] += F12(pp[qq[k]], pp[qq[ko]], kRadius, kSigma);
                  }
                }
              }
            }
          }
        }
      }
    }

    if (fr_s_) {
      for (Size q = 0; q < pp.size(); ++q) {
        Vect r = fr_p_ - pp[q];
        if (r.norm() > kRadius) {
          ff[q] += r * (kPointForceAttractive / std::pow(r.norm(), 3));
        }
      }
    }

    return ff;
  }

 private:
  Vect fr_p_;
  bool fr_s_ = false;
};
