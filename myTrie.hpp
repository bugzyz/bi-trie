#include <stdint.h>
#include <cstddef>
#include <cstdint>
#include <vector>

// for numeric_limits
#include <limits>

// for std::strlen, std::memcpy
#include <cstring>
// for unique_ptr<T>
#include <memory>

// debugging
#include <iostream>
namespace myTrie {
namespace hashRelative {
/**
 * FNV-1a hash
 */
template <class CharT>
struct str_hash {
    std::size_t operator()(const CharT* key, std::size_t key_size) const {
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
    }
};
template <class CharT>
struct str_equal {
    bool operator()(const CharT* key_lhs, std::size_t key_size_lhs,
                    const CharT* key_rhs, std::size_t key_size_rhs) const {
        if (key_size_lhs != key_size_rhs) {
            return false;
        } else {
            return std::memcmp(key_lhs, key_rhs,
                               key_size_lhs * sizeof(CharT)) == 0;
        }
    }
};
}  // namespace hashRelative
}  // namespace myTrie

namespace myTrie {
namespace container {
// todo the hashtable grow policy
/*
template <class CharT, class T,
      class Hash = myTrie::hashRelative::str_hash<CharT>,
      class KeyEqual = myTrie::hashRelative::str_equal<CharT>,
      bool StoreNullTerminator = true, class KeySizeT = std::uint16_t,
      class IndexSizeT = std::uint32_t,
      class GrowthPolicy = tsl::ah::power_of_two_growth_policy<2>>*/
template <class CharT, class T,
          class Hash = myTrie::hashRelative::str_hash<CharT>,
          class KeyEqual = myTrie::hashRelative::str_equal<CharT>,
          bool StoreNullTerminator = true, class KeySizeT = std::uint16_t,
          class IndexSizeT = std::uint32_t>
class array_map {
   public:
    class array_hash_iterator;
    class array_bucket;
    using iterator = array_hash_iterator;

    std::vector<array_bucket> m_buckets;
    std::vector<T> m_values;

    class array_hash_iterator {};
    class array_bucket {
        // store several key like: 5limit\05apple\0
       public:
        CharT* m_buffer;
        class array_bucket_iterator;
        using iterator = array_bucket_iterator;

        // if the length of a key equal to END_OF_BUCKET, it reach the end
        static const KeySizeT END_OF_BUCKET =
            std::numeric_limits<KeySizeT>::max();
        // a key will be followed by one keySizeT of end '\0'
        static const KeySizeT KEY_EXTRA_SIZE = 1;

        iterator end() { return iterator(nullptr); }
        iterator begin() { return iterator(m_buffer); }

        class array_bucket_iterator {
           public:
            array_bucket_iterator(CharT* position) : m_position(position) {}
            CharT* m_position;
        };

        size_t read_key_size(CharT* buffer) {
            // copy the first key_size leading the buffer and return
            KeySizeT key_size;
            std::memcpy(&key_size, buffer, sizeof(key_size));

            return key_size;
        }

        bool is_end_of_bucket(const CharT* buffer) {
            return read_key_size(buffer) == END_OF_BUCKET;
        }

        std::pair<iterator, bool> find_or_end_of_bucket(const CharT* key,
                                                        size_t key_size) {
            if (m_buffer == nullptr) {
                return std::make_pair(end(), false);
            }

            const CharT* buffer_ptr_in_out = m_buffer;
            const bool found =
                find_or_end_of_bucket_impl(key, key_size, buffer_ptr_in_out);

            return std::make_pair(const_iterator(buffer_ptr_in_out), found);
        }

        bool find_or_end_of_bucket_impl(const CharT* key, size_t key_size,
                                        const CharT*& buffer_ptr_in_out) const {
            while (!is_end_of_bucket(buffer_ptr_in_out)) {
                // get key's size
                KeySizeT compared_key_size = read_key_size(buffer_ptr_in_out);
                // skip the key's size
                CharT* buffer_str = buffer_ptr_in_out + sizeof(KeySizeT);
                if (KeyEqual()(buffer_str, compared_key_size, key, key_size)) {
                    return true;
                }

                buffer_ptr_in_out +=
                    entry_size_bytes(buffer_ptr_in_out) / sizeof(CharT);
            }

            return false;
        }

        // todo
        // T is the type of the value
        template <class T>
        static size_t entry_required_bytes(size_t key_size) noexcept {
            return sizeof_in_buff<key_size_type>() +
                   (key_size + KEY_EXTRA_SIZE) * sizeof(CharT) + sizeof(T);
        }
    };

    array_map(size_t bucket_count, float max_load_factor,
              const Hash& hash = Hash()) {
        m_bucket_count = bucket_count;
        m_max_load_factor = max_load_factor;
        m_hash = hash;
    }

    std::pair<iterator, bool> emplace(const CharT* key, size_t key_size, T& v) {
        size_t hashVal = hash_key(key, key_size);
        size_t bucketId = bucket_for_hash(hashVal);
        auto it_find = m_buckets[bucketId].find_or_end_of_bucket(key, key_size);

        // todo: if find

        // todo: bucket grow

        return emplace_impl(bucketId, it_find.first, key, key_size, v);
    }

    std::pair<iterator, bool> emplace_impl(
        std::size_t ibucket, typename array_bucket::iterator end_of_bucket,
        const CharT* key, size_t key_size, T& v) {
        //todo: if reach max

        //todo: if reach capacity

        this->m_values.push_back(v);

        try {
            auto it = m_buckets[ibucket].append(
                end_of_bucket, key, key_size,
                IndexSizeT(this->m_values.size() - 1));
            m_nb_elements++;

            return std::make_pair(
                iterator(m_buckets.begin() + ibucket, it, this), true);
        } catch (...) {
            // Rollback
            this->m_values.pop_back();
            throw;
        }
    }

    size_t bucket_for_hash(std::size_t hash) const {
        return GrowthPolicy::bucket_for_hash(hash);
    }

    std::size_t hash_key(const CharT* key, size_t key_size) const {
        return Hash::operator()(key, key_size);
    }

   public:
    size_t m_bucket_count;
    float m_max_load_factor;
    Hash m_hash;
};
}  // namespace container
}  // namespace myTrie

namespace myTrie {
template <class CharT, class T,
          class Hash = myTrie::hashRelative::str_hash<CharT>,
          class KeySizeT = std::uint16_t>
class htrie_map {
   public:
    class anode;
    class hash_node;
    class trie_node;
    class iterator;

    // htrie_hash
    std::unique_ptr<anode> t_root;
    size_t m_nb_elements;
    float m_max_load_factor;
    size_t m_burst_threshold;
    Hash m_hash;

    // default factor
   public:
    static constexpr float HASH_NODE_DEFAULT_MAX_LOAD_FACTOR = 8.0f;
    static const size_t DEFAULT_BURST_THRESHOLD = 16384;

    static const size_t HASH_NODE_DEFAULT_INIT_BUCKETS_COUNT = 32;
    static const size_t MIN_BURST_THRESHOLD = 4;

   public:
    // htrie_map
    explicit htrie_map(const Hash& hash = Hash())
        : m_hash(hash),
          m_max_load_factor(HASH_NODE_DEFAULT_MAX_LOAD_FACTOR),
          m_burst_threshold(DEFAULT_BURST_THRESHOLD) {
        std::cout << "htrie_map init by the 1 arg consturction " << std::endl;
    }

    T& operator[](const CharT* key) {
        return access_operator(key, std::strlen(key));
    }
    T& operator[](const std::basic_string<CharT>& key) {
        std::cout << "call the [] operator to access" << std::endl;
        return access_operator(key.data(), key.size());
    }

    // htrie_hash
    T& access_operator(const CharT* key, size_t key_size) {
        auto it_find = find(key, key_size);
        if (it_find != end()) {
            return it_find.value();
        } else {
            return insert(key, key_size, T{}).first.value();
        }
    }

    std::pair<iterator, bool> insert(const CharT* key, size_t key_size, T& v) {
        // todo: key size checking

        if (t_root == nullptr) {
            t_root.reset(new hash_node(m_hash, m_max_load_factor));
        }

        return insert_impl(*t_root, key, key_size, v);
    }

    std::pair<iterator, bool> insert_impl(anode& search_start_node,
                                          const CharT* key, size_t key_size,
                                          T& v) {
        anode* current_node = &search_start_node;
        // check whether equal
        for (size_t ikey = 0; ikey < key_size; ikey++) {
            // todo: have content
            if (current_node->isTrieNode()) {
            } else {
                return insert_in_hash_node(current_node->as_hash_node(),
                                           key + ikey, key_size - ikey, v);
            }
        }
    }

    std::pair<iterator, bool> insert_in_hash_node(hash_node& hnode,
                                                  const CharT* key,
                                                  size_t key_size, T& v) {
        // todo burst
        auto it_insert = hnode.array_hash().emplace_ks(key, key_size, v);
        if (it_insert.second) {
            m_nb_elements++;
        }

        return std::make_pair(iterator(hnode, it_insert.first),
                              it_insert.second);
    }

    iterator find(const CharT* key, size_t key_size) {
        if (t_root == nullptr) {
            return end();
        }

        return find_impl(*t_root, key, key_size);
    }

    iterator end() {
        iterator it;
        it.set_as_end_iterator();

        return it;
    }

    // todo
    iterator find_impl(const anode& search_start_node, const CharT* key,
                       size_t key_size) {
        anode* current_node = search_start_node;
    }

    // node relative
    class anode {
       public:
        bool isTrieNode() { return m_node_type == node_type::TRIE_NODE; }

        bool isHashNode() { return m_node_type == node_type::Hash_NODE; }

        enum class node_type : unsigned char { HASH_NODE, TRIE_NODE };
        node_type m_node_type;

        /**
         * If the node has a parent, then it's a descendant of some char.
         *
         * Example:
         *      | ... | a | b | ... | trie_node_1
         *                   \
         *              |...| a | ... | trie_node_2
         *                   /
         *            |array_hash| hash_node_1
         *
         * trie_node_2 is a child of trie_node_1 through 'b', it will have 'b'
         * as m_child_of_char. hash_node_1 is a child of trie_node_2 through
         * 'a', it will have 'a' as m_child_of_char.
         *
         * trie_node_1 has no parent, its m_child_of_char is undeterminated.
         */
        CharT m_child_of_char;
        trie_node* m_parent_node;
    };

    class trie_node : public anode {};

    class hash_node : public anode {
       public:
        typedef myTrie::container::array_map<
            CharT, T, Hash, myTrie::hashRelative::str_equal<CharT>, false,
            KeySizeT, std::uint16_t>
            array_map;
        // todo m_array_hash inside
        hash_node(Hash& hash, float max_load_factor)
            : anode(anode::node_type::hash_node),
              m_array_hash(HASH_NODE_DEFAULT_INIT_BUCKETS_COUNT,
                           HASH_NODE_DEFAULT_MAX_LOAD_FACTOR, hash) {}

        array_map& array_hash() { return m_array_hash; }

        array_map m_array_hash;
    };

    // iterator
    class iterator {
       private:
        trie_node* m_current_trie_node;
        hash_node* m_current_hash_node;

        array_hash_iterator_type m_array_hash_iterator;
        array_hash_iterator_type m_array_hash_end_iterator;

        bool m_read_trie_node_value;
        // TODO can't have void if !IsPrefixIterator, use inheritance
        // typename std::conditional<IsPrefixIterator, std::basic_string<CharT>,
        //                           bool>::type m_prefix_filter;

       public:
        void set_as_end_iterator() {
            m_current_trie_node = nullptr;
            m_current_hash_node = nullptr;
            m_read_trie_node_value = false;
        }
    };
};  // namespace myTrie
}  // namespace myTrie