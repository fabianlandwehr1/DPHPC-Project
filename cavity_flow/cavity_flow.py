# Copyright 2021 ETH Zurich and the NPBench authors. All rights reserved.

import numpy as np


def build_up_b(b, rho, dt, u, v, dx, dy):

    b[1:-1,
      1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                                (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                      ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 *
                      ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                       (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
                      ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    return b


def pressure_poisson(nit, p, dx, dy, b):
    pn = np.empty_like(p)
    pn = p.copy()

    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) - dx**2 * dy**2 /
                         (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])

        p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
        p[-1, :] = 0  # p = 0 at y = 2
    return p


def cavity_flow(nx, ny, nt, nit, u, v, dt, dx, dy, p, rho, nu):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))

    for n in range(nt):
        un = u.copy()
        vn = v.copy()

        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(nit, p, dx, dy, b)

        u[1:-1,
          1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
                   (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                   vn[1:-1, 1:-1] * dt / dy *
                   (un[1:-1, 1:-1] - un[0:-2, 1:-1]) - dt / (2 * rho * dx) *
                   (p[1:-1, 2:] - p[1:-1, 0:-2]) + nu *
                   (dt / dx**2 *
                    (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                    dt / dy**2 *
                    (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))
        # Possible error in reference code:::
        # un used instead of vn in Ln 59, Col 38
        v[1:-1,
          1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx *
                   (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                   vn[1:-1, 1:-1] * dt / dy *
                   (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - dt / (2 * rho * dy) *
                   (p[2:, 1:-1] - p[0:-2, 1:-1]) + nu *
                   (dt / dx**2 *
                    (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                    dt / dy**2 *
                    (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

        u[0, :] = 0
        u[:, 0] = 0
        u[:, -1] = 0
        u[-1, :] = 1  # set velocity on cavity lid equal to 1
        v[0, :] = 0
        v[-1, :] = 0
        v[:, 0] = 0
        v[:, -1] = 0
    return u, v

def initialize(ny, nx):
    u = np.zeros((ny, nx), dtype=np.float64)
    v = np.zeros((ny, nx), dtype=np.float64)
    p = np.zeros((ny, nx), dtype=np.float64)
    dx = 2.0 / (nx - 1)
    dy = 2.0 / (ny - 1)
    dt = .1 / ((nx - 1) * (ny - 1))
    return u, v, p, dx, dy, dt

nx = 4
ny = 4
nit = 2
nt = 4

rho = 1.2
nu = 0.2

u, v, p, dx, dy, dt = initialize(ny,nx)
print(dx, dy, dt)
u, v = cavity_flow(nx, ny, nt, nit, u, v, dt, dx, dy, p, rho, nu)
print(u, v)