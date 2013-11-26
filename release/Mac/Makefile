CC=clang++ -std=c++11 -O3 -DNDEBUG -stdlib=libc++ -I${LIBCXX_ROOT}/include

../../bin/Mac/idock: ../../obj/Mac/utility.o ../../obj/Mac/io_service_pool.o ../../obj/Mac/scoring_function.o ../../obj/Mac/atom.o ../../obj/Mac/receptor.o ../../obj/Mac/ligand.o ../../obj/Mac/random_forest.o ../../obj/Mac/random_forest_x.o ../../obj/Mac/random_forest_.o ../../obj/Mac/log.o ../../obj/Mac/main.o
	$(CC) -o $@ $^ -L${BOOST_ROOT}/lib -lboost_system -lboost_filesystem -lboost_program_options -L${CUDA_ROOT}/lib64 -lcuda -lcurand -lc++abi

../../obj/Mac/main.o: ../../src/main.cpp
	$(CC) -o $@ $< -c -I${BOOST_ROOT} -I${CUDA_ROOT}/include

../../obj/Mac/%.o: ../../src/%.cpp
	$(CC) -o $@ $< -c -I${BOOST_ROOT}

clean:
	rm -f ../../bin/Mac/idock ../../obj/Mac/*.o
