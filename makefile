CUDA := $(shell echo `command -v nvcc`)
MPI := $(shell echo `command -v mpic++`)

all: run_moat moat_a_cpu moat_a_gpu moat_v_serial moat_v_parallel moatsim moatsim_parallel

run_moat: run_moat.cpp
	g++ -Wall -o run_moat run_moat.cpp

moat_a_cpu: mutation_burden_emp.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o moat_a_cpu mutation_burden_emp.cpp variant_permutation_v3.cpp

moat_v_serial: moat_v_pval mutation_burden_emp_v9.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o moat_v_serial mutation_burden_emp_v9.cpp variant_permutation_v3.cpp

moat_a_gpu: mutation_burden_emp.cu variant_permutation_v3.h variant_permutation_v3.cpp
ifndef CUDA
	@echo "No CUDA-capable GPU available for MOAT-a: Can only use CPU version"
else
	nvcc -o moat_a_gpu mutation_burden_emp.cu variant_permutation_v3.cpp
	@echo "CUDA-accelerated MOAT-a is available"
endif

moat_v_parallel: moat_v_pval mutation_burden_emp_v9_mpi.cpp variant_permutation_v3.h variant_permutation_v3.cpp
ifndef MPI
	@echo "No OpenMPI installation detected: Can only use serial version of MOAT-v"
else
	mpic++ -Wall -o moat_v_parallel mutation_burden_emp_v9_mpi.cpp variant_permutation_v3.cpp
	@echo "OpenMPI-accelerated MOAT-v is available"
endif
	
moat_v_pval: p_value_emp.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o p_value_emp p_value_emp.cpp variant_permutation_v3.cpp

moatsim: mutation_burden_emp_v10.cpp variant_permutation_v3.h variant_permutation_v3.cpp
	g++ -Wall -o moatsim mutation_burden_emp_v10.cpp variant_permutation_v3.cpp

moatsim_parallel: mutation_burden_emp_v10_mpi.cpp variant_permutation_v3.h variant_permutation_v3.cpp
ifndef MPI
	@echo "No OpenMPI installation detected: Can only use serial version of MOAT-sim"
else
	mpic++ -Wall -o moatsim_parallel mutation_burden_emp_v10_mpi.cpp variant_permutation_v3.cpp
	@echo "OpenMPI-accelerated MOAT-sim is available"
endif
