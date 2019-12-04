# Using GPUs to solve stiff chemical kinetics problem in reactive flow simulations

The HybridX software is a compressible Navier-Stokes equations solver on structured meshes. This project is to accelerate the software by adding the use of GPUs, which is to be implemented using CUDA.

Team members: Bailin Chen, Jin Dou.



Abstract:

The HybridX software is a compressible Navier-Stokes equations solver on structured meshes. This project is to accelerate the software by adding the use of GPUs, which is to be implemented using CUDA.

HybridX is a well-developed scalable compressible simulation software. It is based on solving the Navier-Stokes equation on structured meshes. Currently HybridX can only utilize the CPUs part of a HPC facility. The goal of this project is to modify the software so that it will be able to take advantage of the GPUs parts of a HPC facility, which is a growing trend nowadays in HPC facility including GUPs powered nodes. 

Due to the nature of the solving fluid dynamics problem, simply moving the calculations from CPUs to GPUs will not accelerate the software because there are too many communication overhead inside the GPUs. To improve the performance overall by using GPUs, the algorithm of the solver is needed to be altered. And due to the scale of the software, the modification of the software will be performed modules by modules of the software. 

For this project in particular, is to move the stiff chemical kinetics problems solving parts to the GPU from the GPUs, as it is proven by recent literature [ref] that it can increase the problem solving performance greatly. 

The schedule of the project is firstly to integrate the GPUs library into the software, which includes adding or changing functions and methods that calls the GPUs library, and changing some of the data structures according so that they are compatible with the GPU library. After that single module testing will be run on the modified code to check its performance and accuracy compared to the original CPU version of the module. If the performance is too poor, which is to be expected, the algorithm of the module will be updated. 
