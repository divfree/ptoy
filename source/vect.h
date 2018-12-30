#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <limits>

template <class T>
T sqr(T a) {
  return a * a;
}

// TODO: move 
template <class Scal>
Scal GetNan() {
  return std::numeric_limits<Scal>::quiet_NaN();
}


template <class Scal, size_t dim_>
class GVect {
 public:
  static constexpr size_t dim = dim_;
  using value_type = Scal;

 private:
  std::array<Scal, dim> aa_;

 public:
  friend void swap(GVect& first, GVect& second) {
    using std::swap;
    swap(first.aa_, second.aa_);
  }
  GVect() = default;
  GVect(const GVect& o) = default;
  size_t size() const {
    return aa_.size();
  }

  // TODO: extend to dim > 2
  explicit GVect(Scal a) : aa_{a, a} {}
  explicit GVect(Scal a0, Scal a1) : aa_({a0, a1}) {}

  template <class, size_t>
  friend class GVect;

  template <class T>
  explicit GVect(const GVect<T, dim>& o) {
    std::copy(o.aa_.begin(), o.aa_.end(), aa_.begin());
  }

  template <class T>
  explicit GVect(const std::vector<T>& v) {
    for (size_t i = 0; i < std::min<size_t>(dim, v.size()); ++i) {
      aa_[i] = static_cast<Scal>(v[i]);
    }
  }
  explicit GVect(const std::array<Scal, dim>& aa) : aa_(aa) {}
  GVect& operator=(GVect other) {
    aa_ = other.aa_;
    return *this;
  }
  Scal& operator[](size_t i) {
    return aa_[i];
  }
  const Scal& operator[](size_t i) const {
    return aa_[i];
  }
  GVect& operator+=(const GVect& vect) {
    for (size_t i = 0; i < dim; ++i) {
      aa_[i] += vect.aa_[i];
    }
    return *this;
  }
  GVect& operator-=(const GVect& other) {
    for (size_t i = 0; i < dim; ++i) {
      aa_[i] -= other.aa_[i];
    }
    return *this;
  }
  GVect& operator*=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      aa_[i] *= k;
    }
    return *this;
  }
  GVect& operator/=(Scal k) {
    for (size_t i = 0; i < dim; ++i) {
      aa_[i] /= k;
    }
    return *this;
  }
  GVect operator+(GVect other) const {
    other += *this;
    return other;
  }
  GVect operator-(GVect other) const {
    GVect tmp(*this);
    tmp -= other;
    return tmp;
  }
  GVect operator-() const {
    GVect tmp(*this);
    for (size_t i = 0; i < dim; ++i) {
      tmp[i] = -tmp[i];
    }
    return tmp;
  }
  GVect operator*(Scal k) const {
    GVect tmp(*this);
    tmp *= k;
    return tmp;
  }
  GVect operator/(Scal k) const {
    GVect tmp(*this);
    tmp /= k;
    return tmp;
  }
  GVect& operator*=(const GVect& other) {
    for (size_t i = 0; i < dim; ++i) {
      aa_[i] *= other.aa_[i];
    }
    return *this;
  }
  GVect operator*(const GVect& other) const {
    GVect tmp(*this);
    tmp *= other;
    return tmp;
  }
  GVect& operator/=(const GVect& other) {
    for (size_t i = 0; i < dim; ++i) {
      aa_[i] /= other.aa_[i];
    }
    return *this;
  }
  GVect operator/(const GVect& other) const {
    GVect tmp(*this);
    tmp /= other;
    return tmp;
  }
  bool operator==(const GVect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(aa_[i] == other.aa_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const GVect& other) const {
    return !(*this == other);
  }
  bool operator<(const GVect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(aa_[i] < other.aa_[i])) {
        return false;
      }
    }
    return true;
  }
  bool operator<=(const GVect& other) const {
    for (size_t i = 0; i < dim; ++i) {
      if (!(aa_[i] <= other.aa_[i])) {
        return false;
      }
    }
    return true;
  }
  bool lexless(const GVect& o) const {
    return aa_ < o.aa_;
  }
  // TODO: remove, replace with GVect(0)
  static const GVect kZero;
  // TODO: remove, replace with GVect(1)
  static GVect GetUnit(size_t i) {
    GVect res = kZero;
    res[i] = 1;
    return res;
  }
  Scal sqrnorm() const {
    Scal res = 0;
    for (size_t i = 0; i < dim; ++i) {
      res += sqr(aa_[i]);
    }
    return res;
  }
  Scal norm() const {
    return std::sqrt(sqrnorm());
  }
  Scal dot(const GVect& other) const {
    Scal sum = 0;
    for (size_t i = 0; i < dim; ++i) {
      sum += aa_[i] * other.aa_[i];
    }
    return sum;
  }
  Scal cross_third(const GVect& other) const {
    return aa_[0] * other.aa_[1] - aa_[1] * other.aa_[0];
  }
  GVect cross(const GVect& other) const {
    const GVect& a = *this;
    const GVect& b = other;
    return GVect(a[1]*b[2]-a[2]*b[1],
                a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0]);
  }
  Scal dist(GVect other) const {
    other -= *this;
    return other.norm();
  }
  Scal sqrdist(GVect o) const {
    o -= *this;
    return o.sqrnorm();
  }
  Scal sum() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r += aa_[i];
    }
    return r;
  }
  Scal prod() const {
    Scal r = aa_[0];
    for (size_t i = 1; i < dim; ++i) {
      r *= aa_[i];
    }
    return r;
  }
  Scal norm1() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r += std::abs(aa_[i]);
    }
    return r;
  }
  Scal norminf() const {
    Scal r = 0;
    for (size_t i = 0; i < dim; ++i) {
      r = std::max(r, std::abs(aa_[i]));
    }
    return r;
  }
  Scal max() const {
    Scal r = aa_[0];
    for (size_t i = 1; i < dim; ++i) {
      r = std::max(r, aa_[i]);
    }
    return r;
  }
  Scal min() const {
    Scal r = aa_[0];
    for (size_t i = 1; i < dim; ++i) {
      r = std::min(r, aa_[i]);
    }
    return r;
  }
  size_t argmax() const {
    size_t r = 0;
    for (size_t i = 1; i < dim; ++i) {
      if (aa_[i] > aa_[r]) {
        r = i;
      }
    }
    return r;
  }
  size_t argmin() const {
    size_t r = 0;
    for (size_t i = 1; i < dim; ++i) {
      if (aa_[i] < aa_[r]) {
        r = i;
      }
    }
    return r;
  }
  GVect abs() const {
    GVect r = *this;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = std::abs(r[i]);
    }
    return r;
  }
  // TODO: revise, may lead to undesired conversion
  template <class T=Scal>
  operator std::vector<T>() const {
    return std::vector<T>(aa_.begin(), aa_.end());
  }
  class LexLess {
   public:
    bool operator()(GVect a, GVect b) const {
      return a.lexless(b);
    }
  };
};

template <class Scal, size_t dim>
const GVect<Scal, dim> GVect<Scal, dim>::kZero =
    GVect<Scal, dim>(static_cast<Scal>(0.));

template <class Scal, size_t dim>
std::ostream& operator<<(std::ostream& out, const GVect<Scal, dim>& vect) {
  out << "(";
  for (size_t i = 0; i < dim; ++i) {
    if (i != 0) {
      out << ",";
    }
    out << vect[i];
  }
  out << ")";
  return out;
}

template <class Scal, size_t dim>
std::istream& operator>>(std::istream& in, GVect<Scal, dim>& vect) {
  for (size_t i = 0; i < dim; ++i) {
    in >> vect[i];
  }
  return in;
}

template <class _Vect>
struct Rect {
  using Vect = _Vect;
  static constexpr size_t dim = Vect::dim;

  Vect lb, rt;
  Rect() {}
  Rect(const Vect& lb, const Vect& rt)
      : lb(lb), rt(rt)
  {}
  bool IsInside(Vect x) const {
    for (size_t i = 0; i < dim; ++i) {
      if (x[i] < lb[i] || rt[i] < x[i]) {
        return false;
      }
    }
    return true;
  }
  Vect GetDimensions() const {
    return rt - lb;
  }
};


