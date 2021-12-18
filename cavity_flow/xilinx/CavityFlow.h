#pragma once

constexpr int NX = 101;
constexpr int NY = 101;

constexpr int NT = 700;
constexpr int NIT = 50;

constexpr float RHO = 1.0;
constexpr float NU = 0.1;

constexpr int D = 4;
constexpr int W = 4;
void CavityFlow(float *u, float *v, float* p);
