#pragma once

#include <cmath>
#include <vector>
#include <cassert>
#include <array>
#include <numeric>

#include "geom.h"


class CellList {
 public:
  // Returns number of cells.
  MIdx GetDims() const {
    return dims_;
  }
  // Returns number of cells.
  Size GetSize() const {
    return dims_.prod();
  }
  // Returns cell index from position
  MIdx GetMIdx(Vect p) const {
    const Scal eps = 1e-5;
    Vect d = (p - dom0_) / (dom1_ - dom0_);
    d[0] = std::max(Scal(0.), std::min(Scal(1.) - eps, d[0]));
    d[1] = std::max(Scal(0.), std::min(Scal(1.) - eps, d[1]));
    return MIdx(d * Vect(dims_));
  }
  Size GetCell(MIdx w) {
    return w[1] * dims_[0] + w[0];
  }
  // pp: positions
  // dims: number of cells
  // dom0,dom1: domain corners
  void Reinit(const Cont<Vect>& pp, MIdx dims, Vect dom0, Vect dom1) {
    dims_ = dims;
    dom0_ = dom0;
    dom1_ = dom1;

    // fill ncc_
    {
      nnc_.clear();
      nnc_.resize(GetSize(), 0);
      for (auto& p : pp) {
        Size c = GetCell(GetMIdx(p));
        ++nnc_[c];
      }
    }

    // fill qq_, kkc_
    {
      kkc_.resize(GetSize());
      // set kkc_ to end of list
      std::partial_sum(nnc_.begin(), nnc_.end(), kkc_.begin());
      assert(kkc_.back() == pp.size());
      nnc_.push_back(pp.size());

      qq_.resize(pp.size());
      ccq_.resize(pp.size());

      // put indices in reversed order
      for (Size q = 0; q < pp.size(); ++q) {
        Size c = GetCell(GetMIdx(pp[q]));
        --kkc_[c];
        qq_[kkc_[c]] = q;
        ccq_[q] = c;
      }
    }
  }
  const Cont<Size>& GetPartIdx() const {
    return qq_;
  }
  const Cont<Size>& GetOffset() const {
    return kkc_;
  }
  const Cont<Size>& GetCellIdx() const {
    return ccq_;
  }

 private:
  MIdx dims_;
  Vect dom0_, dom1_; // domain corners

  Cont<Size> nnc_; // number of particles in cell
  Cont<Size> qq_;  // particle indices
  Cont<Size> kkc_; // offset in qq_ 
  Cont<Size> ccq_; // cell index from particle
};

