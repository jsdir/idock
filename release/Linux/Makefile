CC=g++ -std=c++11 -O3 -DNDEBUG -static

../../bin/Linux/idock: ../../obj/Linux/utility.o ../../obj/Linux/io_service_pool.o ../../obj/Linux/scoring_function.o ../../obj/Linux/atom.o ../../obj/Linux/receptor.o ../../obj/Linux/ligand.o ../../obj/Linux/random_forest.o ../../obj/Linux/random_forest_x.o ../../obj/Linux/random_forest_y.o ../../obj/Linux/log.o ../../obj/Linux/main.o
	$(CC) -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_filesystem -lboost_program_options -L${CUDA_ROOT}/lib64 -lcuda -lcurand

../../obj/Linux/main.o: ../../src/main.cpp
	$(CC) -o $@ $< -c -I${BOOST_ROOT} -I${CUDA_ROOT}/include

../../obj/Linux/%.o: ../../src/%.cpp
	$(CC) -o $@ $< -c -I${BOOST_ROOT}

clean:
	rm -f ../../bin/Linux/idock ../../obj/Linux/*.o
