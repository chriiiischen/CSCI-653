# Using GPUs to solve stiff chemical kinetics problem in reactive flow simulations

The HybridX software is a compressible Navier-Stokes equations solver on structured meshes. This project is to accelerate the software by adding the use of GPUs, which is to be implemented using CUDA.

Team members: Bailin Chen, Jin Dou.

## Abstract:

The topic of our project is Using GPUs to solve stiff chemical kinetics problems in reactive flow simulations, which is actually a part of a larger project: using GPUs to accelerate a software called HybridX. 

HybridX is a well-developed scalable compressible simulation software. It is based on solving the Navier-Stokes equation on structured meshes. Current HybridX can only utilize the CPUs part of a HPC facility. As it is a growing trend nowadays to take advantage of the GPUs parts of a HPC facility, the main task is to move the Solving of Stiff Chemical Kinetics Problems in HybridX from the CPU to the GPU, as it is proven by recent literature [1,2] that it can increase the problem solving performance greatly. 

However, due to the nature of the solving fluid dynamics problem, simply moving the calculations from CPUs to GPUs will not accelerate the software because there are too many communication overhead inside the GPUs. To improve the performance overall by using GPUs, the algorithm of the solver is needed to be altered. And due to the scale of the software, the modification of the software will be performed modules by modules of the software. 

Thus, the project can be generally divided into two stages, the first stage is the integration of CUDA related modules and the second is the modification of HybridX's data structures and algorithms.

For now, the HybridX uses CVODE to solve its ordinary differential equations. CVODE is a solver for stiff and nonstiff ordinary differential equation (ODE) systems (initial value problem) given in explicit form y’ = f(t,y). NVECTOR is the Key NVECTOR_CUDA module is an experimental NVECTOR implementation in the cuda language. It works as an adapter between CVODE and CUDA. So for the first stage, we need to integrate the NVECTOR_CUDA module into the software, which includes adding or changing functions and methods to let the computation part can run seccuessfully in GPUS.

For the second stage, to reduce the communication overhead, we firstly have to change some of the data structures according so that they can be kept on the GPU. After that single module testing will be run on the modified code to check its performance and accuracy compared to the original CPU version of the module. If the performance is too poor, which is to be expected, the algorithm of the module will be updated. 

## Reference:

[1] Stone C P, Alferman A T, Niemeyer K E. Accelerating finite-rate chemical kinetics with coprocessors: comparing vectorization methods on GPUs, MICs, and CPUs[J]. Computer Physics Communications, 2018, 226: 18-29.

[2] Niemeyer K E, Sung C J. Accelerating moderately stiff chemical kinetics in reactive-flow simulations using GPUs[J]. Journal of Computational Physics, 2014, 256: 854-871.
