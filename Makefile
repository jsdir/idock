CC=clang++ -std=c++11 -O2
NVCC=nvcc -Xptxas=-v -use_fast_math

all: bin/idock bin/idock.fatbin

bin/idock: obj/utility.o obj/io_service_pool.o obj/scoring_function.o obj/atom.o obj/receptor.o obj/ligand.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/log.o obj/main.o
	$(CC) -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_program_options -lboost_filesystem -L${CUDA_ROOT}/lib64 -lcuda -lcurand

obj/main.o: src/main.cpp
	$(CC) -o $@ $< -c -I${BOOST_ROOT} -I${CUDA_ROOT}/include

obj/%.o: src/%.cpp
	$(CC) -o $@ $< -c -I${BOOST_ROOT}

bin/%.fatbin: src/%.cu
	$(NVCC) -o $@ $< -fatbin -gencode arch=compute_11,code=compute_11

clean:
	rm -f bin/idock bin/idock.fatbin obj/*.o
