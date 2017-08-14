CUDA := $(shell echo `command -v nvcc`)
MPI := $(shell echo `command -v mpic++`)

all: run_moat moat_a_cpu moat_a_gpu moat_v_serial moat_v_parallel moatsim moatsim_parallel

run_moat: run_moat.cpp
	g++ -Wall -o run_moat run_moat.cpp

# moat_a_cpu: moat_a.cpp variant_permutation_v3.h variant_permutation_v3.cpp
# 	g++ -Wall -o moat_a_cpu moat_a.cpp variant_permutation_v3.cpp

moat_v_serial: moat_v_pval moat_v.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o moat_v_serial moat_v.cpp variant_permutation_v3.cpp

moat_a_gpu: moat_a.cu variant_permutation_v3.h variant_permutation_v3.cpp
ifndef CUDA
	@echo "No CUDA-capable GPU available for MOAT-a: Can only use CPU version"
else
	nvcc -o moat_a_gpu moat_a.cu variant_permutation_v3.cpp
	nvcc -o moat_a_v2 moat_a_v2.cu variant_permutation_v3.cpp
	@echo "CUDA-accelerated MOAT-a is available"
endif

moat_a_v2: moat_a_v2.cu variant_permutation_v3.h variant_permutation_v3.cpp
ifndef CUDA
	@echo "No CUDA-capable GPU available for MOAT-a v2: Can only use CPU version"
else
	nvcc -o moat_a_v2 moat_a_v2.cu variant_permutation_v3.cpp
	@echo "CUDA-accelerated MOAT-a v2 is available"
endif

moat_v_parallel: moat_v_pval moat_v_mpi.cpp variant_permutation_v3.h variant_permutation_v3.cpp
ifndef MPI
	@echo "No OpenMPI installation detected: Can only use serial version of MOAT-v"
else
	mpic++ -Wall -o moat_v_parallel moat_v_mpi.cpp variant_permutation_v3.cpp
	@echo "OpenMPI-accelerated MOAT-v is available"
endif
	
moat_v_pval: p_value_emp.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o p_value_emp p_value_emp.cpp variant_permutation_v3.cpp

moatsim: moatsim.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o moatsim moatsim.cpp variant_permutation_v3.cpp

moatsim_parallel: moatsim_mpi.cpp variant_permutation_v3.h variant_permutation_v3.cpp
ifndef MPI
	@echo "No OpenMPI installation detected: Can only use serial version of MOAT-sim"
else
	mpic++ -Wall -o moatsim_parallel moatsim_mpi.cpp variant_permutation_v3.cpp
	@echo "OpenMPI-accelerated MOATsim is available"
endif
