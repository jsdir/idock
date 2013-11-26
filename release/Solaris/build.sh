rm -f ../../bin/Solaris/idock ../../obj/Solaris/*.o
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/utility.o ../../src/utility.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/io_service_pool.o ../../src/io_service_pool.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/scoring_function.o ../../src/scoring_function.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/atom.o ../../src/atom.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/receptor.o ../../src/receptor.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/ligand.o ../../src/ligand.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/grid_map_task.o ../../src/grid_map_task.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/monte_carlo_task.o ../../src/monte_carlo_task.cpp -c -I${BOOST_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../obj/Solaris/main.o ../../src/main.cpp -c -I${BOOST_ROOT} -I${CUDA_ROOT}
g++ -std=c++0x -O3 -DNDEBUG -o ../../bin/Solaris/idock ../../obj/Solaris/utility.o ../../obj/Solaris/io_service_pool.o ../../obj/Solaris/scoring_function.o ../../obj/Solaris/atom.o ../../obj/Solaris/receptor.o ../../obj/Solaris/ligand.o ../../obj/Solaris/random_forest.o ../../obj/Solaris/random_forest_x.o ../../obj/Solaris/random_forest_y.o ../../obj/Solaris/log.o ../../obj/Solaris/main.o -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_filesystem -lboost_program_options -L${CUDA_ROOT} -lcuda -lcurand
