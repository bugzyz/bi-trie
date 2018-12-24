// decide test which version of trie-map
// #define TEST_CUCKOOHASH
// #define TEST_GROWCUCKOOHASH
// #define TEST_HAT
#define SHRINK_TEST_GROWCUCKOOHASH

// decide the growing cuckoo hash is rehash first or expand first
// #define REHASH_BEFORE_EXPAND

// decide grow the associativity or bucket
// #define GROW_BUCKET

#ifdef TEST_CUCKOOHASH
#include "impl_cuckoo_hash/cuckoo_trie.hpp"
#endif
#ifdef TEST_HAT
#include "impl_hat_trie/hat_trie.hpp"
#endif
#ifdef TEST_GROWCUCKOOHASH
#ifndef GROW_BUCKET
#include "impl_growing_cuckoo_hash/grow_associativity/grow_cuckoo_trie_ass.hpp"
#else
#include "impl_growing_cuckoo_hash/grow_bucket/grow_cuckoo_trie_buc.hpp"
#endif
#endif

#ifdef SHRINK_TEST_GROWCUCKOOHASH
#include "impl_growing_cuckoo_hash_shrink/grow_cuckoo_trie_ass_shrnk.hpp"
#endif
