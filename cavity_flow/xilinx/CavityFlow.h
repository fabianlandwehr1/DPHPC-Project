#pragma once

#include "hlslib/xilinx/DataPack.h"

constexpr int NY = 32;
constexpr int NX = 32;

constexpr int NT = 25;
constexpr int NIT = 5;

constexpr float RHO = 1.0;
constexpr float NU = 0.1;

constexpr int D = 4;
constexpr int W = 4;

using VecW = hlslib::DataPack<float, W>;
void CavityFlow(float *u, float *v, float* p);
