#pragma once

constexpr int NX = 4;
constexpr int NY = 4;

constexpr int NT = 4;
constexpr int NIT = 2;

constexpr float RHO = 1.2;
constexpr float NU = 0.2;

constexpr int D = 4;
constexpr int W = 4;
void CavityFlow(float *u, float *v, float* p);
