rm -f ../bin/Solaris/idock ../obj/Solaris/*.o
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/utility.o ../src/utility.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/io_service_pool.o ../src/io_service_pool.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/scoring_function.o ../src/scoring_function.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/atom.o ../src/atom.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/receptor.o ../src/receptor.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/ligand.o ../src/ligand.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/grid_map_task.o ../src/grid_map_task.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/monte_carlo_task.o ../src/monte_carlo_task.cpp -I${BOOST_ROOT} -c
g++ -std=c++0x -O3 -DNDEBUG -pthreads -o ../obj/Solaris/main.o ../src/main.cpp -I${BOOST_ROOT} -c
g++ -O3 -DNDEBUG -o ../bin/Solaris/idock ../obj/Solaris/utility.o ../obj/Solaris/io_service_pool.o ../obj/Solaris/scoring_function.o ../obj/Solaris/atom.o ../obj/Solaris/receptor.o ../obj/Solaris/ligand.o ../obj/Solaris/random_forest.o ../obj/Solaris/random_forest_x.o ../obj/Solaris/random_forest_y.o ../obj/Solaris/log.o ../obj/Solaris/main.o -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_filesystem -lboost_program_options -L${CUDA_ROOT} -lcuda -lcurand
