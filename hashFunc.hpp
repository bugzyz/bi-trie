// size_t
#include <cstdint>

// for std::memcmp, std::strlen
#include <cstring>

#include <map>
#include <vector>

#include <string>

// numeric_limit
#include <limits>

namespace myTrie {
namespace hashRelative {

// For fasthash64:
static inline uint64_t mix(uint64_t h) {
    h ^= h >> 23;
    h *= 0x2127599bf4325c37ULL;
    h ^= h >> 47;
    return h;
}

// A default hash function:
uint64_t fasthash64(const char* buf, size_t len, uint64_t seed) {
    uint64_t const m = 0x880355f21e6d1965ULL;
    uint64_t const* pos = (uint64_t const*)buf;
    uint64_t const* end = pos + (len / 8);
    const unsigned char* pos2;
    uint64_t h = seed ^ (len * m);
    uint64_t v;

    while (pos != end) {
        v = *pos++;
        h ^= mix(v);
        h *= m;
    }

    pos2 = (const unsigned char*)pos;
    v = 0;

    switch (len & 7) {
        case 7:
            v ^= (uint64_t)pos2[6] << 48;
        case 6:
            v ^= (uint64_t)pos2[5] << 40;
        case 5:
            v ^= (uint64_t)pos2[4] << 32;
        case 4:
            v ^= (uint64_t)pos2[3] << 24;
        case 3:
            v ^= (uint64_t)pos2[2] << 16;
        case 2:
            v ^= (uint64_t)pos2[1] << 8;
        case 1:
            v ^= (uint64_t)pos2[0];
            h ^= mix(v);
            h *= m;
    }

    return mix(h);
}

template <class CharT>
std::size_t hash(const CharT* key, std::size_t key_size, size_t hashType) {
    uint64_t hash;
    if (hashType == 1) {
        hash = fasthash64(key, key_size, 0xdeadbeefdeadbeefULL);
        return hash;
    } else {
        hash = fasthash64(key, key_size, 0xabcdefabcdef1234ULL);
        return hash;
    }
}

template <class CharT>
bool keyEqual(const CharT* key_lhs, std::size_t key_size_lhs,
              const CharT* key_rhs, std::size_t key_size_rhs) {
    if (key_size_lhs == 0 && key_size_rhs == 0) {
        return true;
    }
    if (key_size_lhs != key_size_rhs) {
        return false;
    } else {
        return std::memcmp(key_lhs, key_rhs, key_size_lhs * sizeof(CharT)) == 0;
    }
}
}  // namespace hashRelative
}  // namespace myTrie