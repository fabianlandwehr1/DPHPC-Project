# Hand-optimized FPGA Implementations of NPBench Kernels

_Field Programmable Gate Arrays (FPGA) offer an alternate computation platform in addition to traditional CPU- and GPU-based systems. Unfortunately, it is difficult to fully realize the advantages of FPGAs due to programming challenges characterized by, for example, the lack of automatic memory management and programmer-transparent caches. To this end, various frameworks aim to compile code written in a high-level language such as Python to an FPGA-compatible bitstream. Evaluating these frameworks involves measuring how close to optimal the performance of an automatically generated implementation is, introducing the need for hand-optimized code as reference. This project provides such FPGA implementations of five numerical kernels taken from the [NPBench](https://github.com/spcl/npbench) benchmark suite. We experimentally evaluate our implementations using synthesis, emulation and hardware execution on a state-of-the-art Xilinx FPGA. Our implementations achieve up to 10.7x higher performance than comparable CPU- and GPU-based versions._

## Benchmarks

- **Azimint:** [baseline](https://github.com/fabianlandwehr1/DPHPC-Project/blob/azimint/base/azimint/xilinx/Azimint.cpp) - [opt1](https://github.com/fabianlandwehr1/DPHPC-Project/blob/azimint/opt1/azimint/xilinx/Azimint.cpp) - [opt2](https://github.com/fabianlandwehr1/DPHPC-Project/blob/azimint/opt2/azimint/xilinx/Azimint.cpp)
- **Durbin:** [baseline](https://github.com/fabianlandwehr1/DPHPC-Project/blob/durbin/base/durbin/xilinx/Durbin.cpp) - [opt1](https://github.com/fabianlandwehr1/DPHPC-Project/blob/durbin/opt1/durbin/xilinx/Durbin.cpp) - [opt2](https://github.com/fabianlandwehr1/DPHPC-Project/blob/durbin/opt2/durbin/xilinx/Durbin.cpp)
- **Gram-Schmidt:** [baseline](https://github.com/fabianlandwehr1/DPHPC-Project/blob/main/gramschmidt/xilinx/naive/GramSchmidt1.cpp) - [opt1](https://github.com/fabianlandwehr1/DPHPC-Project/blob/main/gramschmidt/xilinx/third_implementation/GramSchmidt3.cpp)
- **Cavity Flow:** [baseline](https://github.com/fabianlandwehr1/DPHPC-Project/blob/cavityflowNaive/cavity_flow/xilinx/CavityFlow.cpp) - [opt1](https://github.com/fabianlandwehr1/DPHPC-Project/blob/cavityFlowPipelined/cavity_flow/xilinx/CavityFlow.cpp)
- **Conv2D:** [baseline](https://github.com/fabianlandwehr1/DPHPC-Project/blob/main/conv2d_bias/xilinx/naive/conv2d1.cpp) - [opt1](https://github.com/fabianlandwehr1/DPHPC-Project/blob/main/conv2d_bias/xilinx/buffered_and_vectorized/conv2d2.cpp) - [opt2](https://github.com/fabianlandwehr1/DPHPC-Project/blob/main/conv2d_bias/xilinx/with_loop_reordering_and_fully_pipelined/conv2d3.cpp)

## Building

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
make TestAzimintXilinx
azimint/TestAzimintXilinx
```

Synthesize HLS for Xilinx:

```shell
make synthesize_azimint
```
