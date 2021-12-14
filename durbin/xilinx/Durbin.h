#pragma once

#include "hlslib/xilinx/DataPack.h"


constexpr int N = 1024;
constexpr int num_streams = 128;
void DurbinMod(float const r[], float y[]);
void Durbin(float const r[], float y[]);

constexpr int D = N /num_streams;
constexpr int W = 4;
using Vec_t = hlslib::DataPack<float, W>;
