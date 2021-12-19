#pragma once

constexpr int NX = 121;
constexpr int NY = 121;

constexpr int NT = 50;
constexpr int NIT = 10;

constexpr float RHO = 1.0;
constexpr float NU = 0.1;

constexpr int D = 4;
constexpr int W = 4;
void CavityFlow(float *u, float *v, float* p);
