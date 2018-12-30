#pragma once

#include <cmath>
#include <vector>
#include <cassert>
#include <array>

#include "geom.h"


class CellList {
 public:
  // Returns number of cells.
  MIdx GetSize() const {
    return size_;
  }
  // Returns cell index from position
  MIdx PosToCell(Vect p) const {
    const Scal eps = 1e-5;
    Vect d = (p - dom0_) / (dom1_ - dom0_);
    d[0] = std::max(Scal(0.), std::min(Scal(1.) - eps, d[0]));
    d[1] = std::max(Scal(0.), std::min(Scal(1.) - eps, d[1]));
    return MIdx(d * Vect(size_));
  }
  // pp: positions
  // size: number of cells
  // dom0,dom1: domain corners
  void Reinit(const Cont<Vect>& pp, MIdx size, Vect dom0, Vect dom1) {
    size_ = size;
    dom0_ = dom0;
    dom1_ = dom1;
  }

 private:
  MIdx size_;
  Vect dom0_, dom1_; // domain corners

  Cont<Size> cq_; // cell to particle
  Cont<Size> qc_; // particle to cell
};

