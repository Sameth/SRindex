benchmarks: graph.cpp sr-index.cpp benchmark/query.cpp
	g++ -Wall -Wextra -O2 -DNDEBUG -I ~/include -L ~/lib -g -std=c++11 -o benchmark/query.bin sr-index.cpp graph.cpp benchmark/query.cpp -lsdsl -ldivsufsort -ldivsufsort64 -lboost_system -lboost_filesystem
	g++ -Wall -Wextra -O2 -DNDEBUG -I ~/include -L ~/lib -g -std=c++11 -o benchmark/test.bin sr-index.cpp graph.cpp benchmark/test.cpp -lsdsl -ldivsufsort -ldivsufsort64 -lboost_system -lboost_filesystem
