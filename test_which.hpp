// #define TEST_CUCKOOHASH
// #define TEST_GROWCUCKOOHASH
#define TEST_HAT

#ifdef TEST_CUCKOOHASH
#include "impl_cuckoo_hash/myTrie.hpp"
#else
#ifdef TEST_HAT
#include "impl_hat_trie/myTrie.hpp"
#else
#include "impl_growing_cuckoo_hash/myTrie.hpp"
#endif
#endif