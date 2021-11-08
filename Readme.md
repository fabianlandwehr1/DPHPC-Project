# Repository for the DPHPC Project

Clone the repository (including definelicht/hlslib submodule):

```shell
git clone --recursive git@github.com:fabianlandwehr1/DPHPC-Project.git
```

Configure build system using CMake:

```shell
mkdir build
cd buid
cmake ..
```

Build and run test (example):

```shell
make TestCholeskyXilinx
cholesky/TestCholeskyXilinx
```

Synthesize HLS for Xilinx:

```shell
make synthesize_cholesky
```
