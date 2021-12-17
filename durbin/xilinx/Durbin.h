#pragma once

#include "hlslib/xilinx/DataPack.h"


constexpr int N = 1024;
constexpr int num_streams = 32;
constexpr int D = N /num_streams;
constexpr int W = 4;
constexpr int numItersMain = N/D + 1;
constexpr int numItersFirstIter = N/D;

using Vec_t = hlslib::DataPack<float, W>;

void DurbinMod(float const r[], float y[]);
void Durbin(float const r[], float y[]);

