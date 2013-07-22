CC = clang++ -O3 -DNDEBUG -std=c++11
NVCC = nvcc -arch=sm_11 -Xptxas=-v -ftz=true -prec-div=false -prec-sqrt=false -use_fast_math# -maxrregcount=N

bin/idock: obj/utility.o obj/thread_pool.o obj/scoring_function.o obj/atom.o obj/receptor.o obj/ligand.o obj/kernel.o obj/main.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_program_options -lboost_filesystem -L /opt/cuda/lib64 -lcudart -lcurand

obj/kernel.o: src/kernel.cu
	$(NVCC) -o $@ $< -c -I${CUDA_ROOT}/samples/common/inc

obj/%.o: src/%.cpp 
	$(CC) -o $@ $< -c

clean:
	rm -f bin/idock obj/*.o
