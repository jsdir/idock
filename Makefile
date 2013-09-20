CC=clang++ -std=c++11 -O3
NVCC=nvcc -gencode arch=compute_11,code=sm_11 -Xptxas=-v -ftz=true -prec-div=false -prec-sqrt=false -use_fast_math# -maxrregcount=N -g -G -Xcompiler -rdynamic -lineinfo -src-in-ptx -cubin -gencode arch=compute_10,code=compute_10

bin/idock: obj/utility.o obj/thread_pool.o obj/scoring_function.o obj/atom.o obj/receptor.o obj/ligand.o obj/cu_engine.o obj/cu_mc_kernel.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/main.o #obj/cl_engine.o obj/cl_mc_kernel.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_program_options -lboost_filesystem -L${CUDA_ROOT}/lib64 -lcudart -lcurand

obj/%.o: src/%.cu
	$(NVCC) -o $@ $< -c -I${CUDA_ROOT}/samples/common/inc

obj/cu_%.o: src/cu_%.cpp
	$(CC) -o $@ $< -c -I${CUDA_ROOT}/include -I${CUDA_ROOT}/samples/common/inc

obj/%.o: src/%.cpp
	$(CC) -o $@ $< -c

clean:
	rm -f bin/idock obj/*.o
