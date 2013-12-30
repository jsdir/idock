CC=clang++ -std=c++11 -O2

all: combinelog combinelog2 extractelitists extractmodel parsetime pdbqt2csv rmsd

combinelog: combinelog.cpp
	$(CC) -o $@ $< -lboost_system -lboost_filesystem

combinelog2: combinelog2.cpp
	$(CC) -o $@ $< -lboost_system -lboost_filesystem

extractelitists: extractelitists.cpp
	$(CC) -o $@ $< -lboost_system -lboost_filesystem

extractmodel: extractmodel.cpp
	$(CC) -o $@ $< -lboost_system -lboost_filesystem

parsetime: parsetime.cpp
	$(CC) -o $@ $< -lboost_system -lboost_filesystem

pdbqt2csv: pdbqt2csv.cpp
	$(CC) -o $@ $< -lboost_system -lboost_filesystem

rmsd: rmsd.cpp
	$(CC) -o $@ $<
