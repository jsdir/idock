CC=clang++ -std=c++11 -O2

bin/idock: obj/scoring_function.o obj/box.o obj/quaternion.o obj/io_service_pool.o obj/safe_counter.o obj/atom.o obj/receptor.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/log.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/main.o
	$(CC) -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_program_options -lboost_filesystem

obj/%.o: src/%.cpp 
	$(CC) -o $@ $< -c -I${BOOST_ROOT}

clean:
	rm -f bin/idock obj/*.o
