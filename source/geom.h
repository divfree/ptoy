#pragma once

#include <vector>

#include "vect.h"
#include "aligned_allocator.hpp"

using Scal = float;
using Vect = GVect<Scal, 2>;
template <class T>
using Cont = std::vector<T, hpc15::aligned_allocator<T,64>>;
using Size = unsigned short;
using MIdx = GVect<Size, 2>;


const Scal kRadius = 0.02;
const Scal kSigma = 1.;
const Scal kSigmaWall = 1.;
const Scal kSigmaBond = 1e5;
const Scal kSigmaPick = 1e3;
const Scal kSigmaPortalEdge = 1;
const Scal kRadiusPortalEdge = 3. * kRadius;
const Scal kMass = kRadius * kRadius * 100.;
const Scal kPointForce = 0.1;
const Scal kPointForceAttractive = 0.1;
const Scal kDissipation = .01;
const Scal kTimeStep = 0.0003;
const Scal kBlockSize = 4. * kRadius;
const Scal kGravity = 10.;
const Scal kPortalThickness = 0.02;
const Scal kVelocityLimit = 10.;
