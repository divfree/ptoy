#pragma once

#include "vect.h"
#include "aligned_allocator.hpp"

using Scal = float;
using Vect = GVect<Scal, 2>;
template <class T>
using Cont = std::vector<T, hpc15::aligned_allocator<T,64>>;
using Size = unsigned short;
using MIdx = GVect<Size, 2>;
