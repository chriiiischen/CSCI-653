# Using GPUs to solve stiff chemical kinetics problem in reactive flow simulations

### Introduction
The HybridX software is a compressible Navier-Stokes equations solver on structured meshes. This project is to accelerate the software by adding the use of GPUs, which is to be implemented using CUDA.

### Team members 
Bailin Chen, Jin Dou.

## Abstract

The topic of our project is Using GPUs to solve stiff chemical kinetics problems in reactive flow simulations, which is actually a part of a larger project: using GPUs to accelerate a software called HybridX. 

HybridX is a well-developed scalable compressible simulation software. It is based on solving the Navier-Stokes equation on structured meshes. Current HybridX can only utilize the CPUs part of a HPC facility. As it is a growing trend nowadays to take advantage of the GPUs parts of a HPC facility, the main task is to move the Solving of Stiff Chemical Kinetics Problems in HybridX from the CPU to the GPU, as it is proven by recent literature [1,2] that it can increase the problem solving performance greatly. 

However, due to the nature of the solving fluid dynamics problem, simply moving the calculations from CPUs to GPUs will not accelerate the software because there are too many communication overhead inside the GPUs. To improve the performance overall by using GPUs, the algorithm of the solver is needed to be altered. And due to the scale of the software, the modification of the software will be performed modules by modules of the software. 

Thus, the project can be generally divided into two stages, the first stage is the integration of CUDA related modules and the second is the modification of HybridX's data structures and algorithms.

For now, the HybridX uses CVODE to solve its ordinary differential equations. CVODE is a solver for stiff and nonstiff ordinary differential equation (ODE) systems (initial value problem) given in explicit form y’ = f(t,y). Because the sundials solvers of CVODE are written in a data-independent manner, all of them are operated on generic vectors called NVECTOR, which is in charge of key data storage. To make NVECTOR can also stay on the CUDA based GPU. CVODE provides a experimental NVECTOR implementation in the cuda language called NVECTOR_CUDA. It works as an adapter between CVODE and CUDA. So for the first stage, we need to integrate the NVECTOR_CUDA module into the software, which includes adding or changing functions and methods to let the computation part can run seccuessfully in GPUS.

For the second stage, to reduce the communication overhead, we firstly have to change some of the data structures according so that they can be kept on the GPU. After that single module testing will be run on the modified code to check its performance and accuracy compared to the original CPU version of the module. If the performance is too poor, which is to be expected, the algorithm of the module will be updated. 


### Up to now

Up to now, the CUDA version of NVECTOR was integrated successfully with the target component. The main file that was modified was the Integrator of the component. Both the original file (mathCVodesIntegrator.cpp) and the modified file (mathCVodesIntegrator_GPU.cu) are provided in this repository. The Integrator is a hidden layer which is called by the “cvode_gas” class in the programme. So when running the test, the test file (testphysics_cvode_gas_stand_alone.cpp, provided in the repository) does not need to be modified (expect changing the extension to .cu). The CUDA version of NVECTOR provides the options to use managed or unmanaged CUDA memory. In here, the unmanaged option was chosen so that we can have higher control of the programme. By doing so, optimizing the performance in the next step will be easier.

#### Current Testing and Experiment result:
The testing (testphysics_cvode_gas_stand_alone.cpp) that was used is a stand alone test. This means that the programme only tests the target component. It uses a single cell (single integration equation) with a small number of chemistry species and first initialize/reinitialize, then integrate in each step for 100 steps. The purpose of this preliminary test is to validate the accuracy and the consistency of the modification. The results indicated that the modified component gave the exact same results as the original component gave. However, in terms of the performance, the cpu version performs better than the gpu version. The time spent on each step for each version was recorded and shown as follows:

Single step time:
CPU: 4.02e-06   GPU:1.54e-04

This result was expected due to the fact that the optimization step has not been performed. Furthermore, the calculation in the current test is too light to see the speedup using GPUs. 

### Current work:
1. Add loop in the test programme to build a 2-D grids of cells where each cell host an integration process and change the profile with larger number of chemistry species to increase the calculation intensity. The reason to do this is to validate the benefit of using GPUs in the calculation.
2. Apply parallelism to the loop (may need to modify the data structure). The reason to do this is to increase the performance of the programme. 

### Next Step:
1. Do more experiments based on the improved test programme.
2. Further modify the mathCVodesIntegrator_GPU.cu component based on the experiment result.


## Reference

[1] Stone C P, Alferman A T, Niemeyer K E. Accelerating finite-rate chemical kinetics with coprocessors: comparing vectorization methods on GPUs, MICs, and CPUs[J]. Computer Physics Communications, 2018, 226: 18-29.

[2] Niemeyer K E, Sung C J. Accelerating moderately stiff chemical kinetics in reactive-flow simulations using GPUs[J]. Journal of Computational Physics, 2014, 256: 854-871.
