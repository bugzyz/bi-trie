// decide test which version of trie-map
#define TEST_CUCKOOHASH
// #define TEST_GROWCUCKOOHASH
// #define TEST_HAT

// decide the growing cuckoo hash is rehash first or expand first
#define REHASH_BEFORE_EXPAND

#ifdef TEST_CUCKOOHASH
#include "impl_cuckoo_hash/cuckoo_trie.hpp"
#endif
#ifdef TEST_HAT
#include "impl_hat_trie/hat_trie.hpp"
#endif
#ifdef TEST_GROWCUCKOOHASH
#include "impl_growing_cuckoo_hash/grow_cuckoo_trie.hpp"
#endif