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

template <class CharT>
std::size_t hash(const CharT* key, std::size_t key_size) {
    static const std::size_t init = std::size_t(
        (sizeof(std::size_t) == 8) ? 0xcbf29ce484222325 : 0x811c9dc5);
    static const std::size_t multiplier =
        std::size_t((sizeof(std::size_t) == 8) ? 0x100000001b3 : 0x1000193);

    std::size_t hash = init;
    for (std::size_t i = 0; i < key_size; ++i) {
        hash ^= key[i];
        hash *= multiplier;
    }

    return hash;

    //todo: test if the hash func should be changed
    // std::string s(key, key_size);
    // return std::hash<std::string>{}(s);
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