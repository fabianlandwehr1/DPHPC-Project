#pragma once

constexpr int N = 1000000;
constexpr int npt = 1000;
constexpr int d = 10;
constexpr int D = npt / d;
constexpr int W = 4;

void Azimint(float const *data, float const *radius, float *res);
