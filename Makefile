test: 
	g++ gtest/test_bi_trie.cc -std=c++17 -g -O0 -I boost_1_67_0 -I googletest/googletest/include -L googletest/lib -pthread -lgtest -lgtest_main

batch_test:
	g++ test_file/batch_test.cpp -std=c++17 -g -O0 -I boost_1_67_0
	