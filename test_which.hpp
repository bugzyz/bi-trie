// decide test which version of trie-map
// #define TEST_HAT
#define TEST_HASH
// #define TEST_CUCKOOHASH
// #define TEST_GROW_CUCKOOHASH_ASS
// #define TEST_GROW_CUCKOOHASH_BUC

// decide the growing cuckoo hash is rehash first or expand first
// #define REHASH_BEFORE_EXPAND

// decide whether test yago dataset
// #define TEST_YAGO

#ifdef TEST_HAT
#include "unified_impl/1_tessil_hat_impl.hpp"
#endif
#ifdef TEST_HASH
#include "unified_impl/2_prototype_shrink.hpp"
#endif
#ifdef TEST_CUCKOOHASH
#include "unified_impl/3_prototype_cuckoo_shrink.hpp"
#endif
#ifdef TEST_GROW_CUCKOOHASH_ASS
#include "unified_impl/4_prototype_cuckoo_grow_ass_shrink.hpp"
#endif
#ifdef TEST_GROW_CUCKOOHASH_BUC
#include "unified_impl/5_prototype_cuckoo_grow_buc_shrink.hpp"
#endif
