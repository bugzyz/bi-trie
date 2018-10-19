#include <cstdint>
// for std::strlen
#include <cstring>
#include <map>
#include <vector>

#include <string>

// numeric_limit
#include <limits>

// debug
#include <iostream>

#define BUCKET_INIT_COUNT 5

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

    return hash % BUCKET_INIT_COUNT;
}

template <class CharT>
bool keyEqual(const CharT* key_lhs, std::size_t key_size_lhs,
              const CharT* key_rhs, std::size_t key_size_rhs) {
    if (key_size_lhs != key_size_rhs) {
        return false;
    } else {
        return std::memcmp(key_lhs, key_rhs, key_size_lhs * sizeof(CharT)) == 0;
    }
}
}  // namespace hashRelative
}  // namespace myTrie

namespace myTrie {
// charT = char, T = value type, keysizeT = the type describe keysize
template <class CharT, class T, class KeySizeT = std::uint16_t>
class htrie_map {
   public:
    // TODO: burst_threshold should be set to be greater than 26 for 26 alaphbet
    // and ", @ .etc.
    static const size_t DEFAULT_BURST_THRESHOLD = 2;
    size_t burst_threshold;

    enum class node_type : unsigned char { HASH_NODE, TRIE_NODE };
    class trie_node;
    class hash_node;
    class array_bucket;

    class anode {
       public:
        node_type _node_type;
        trie_node* parent;
        CharT myChar;

        bool isHashNode() { return _node_type == node_type::HASH_NODE; }
        bool isTrieNode() { return _node_type == node_type::TRIE_NODE; }
        void setParent(trie_node* p) { parent = p; }
        void setParent(CharT c) { myChar = c; }
    };
    class trie_node : public anode {
       public:
        // store the suffix of hash_node or trie_node
        std::map<CharT, anode*> childs;

        trie_node(CharT c) {
            anode::_node_type = node_type::TRIE_NODE;
            anode::myChar = c;
        }

        std::map<CharT, anode*> getChildsMap() { return childs; }

        anode* findChildNode(CharT c) {
            auto found = childs.find(c);
            if (found != childs.end()) {
                return found->second;
            } else {
                return nullptr;
            }
        }

        void addChildNode(CharT c, anode* node) { childs[c] = node; }

        void addChileNode(std::vector<std::pair<CharT, anode*>> newChilds) {
            for (int i = 0; i != newChilds.size(); i++) {
                childs[newChilds[i].first] = newChilds[i].second;
            }
        }
    };

    class hash_node : public anode {
       public:
        std::vector<array_bucket> kvs;
        size_t burst_threshold;

        hash_node(size_t customized_burst_threshold) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = nullptr;
            kvs = std::vector<array_bucket>(BUCKET_INIT_COUNT);
            burst_threshold = customized_burst_threshold;
        }

        bool need_burst() {
            size_t node_element_count;
            for (auto it = kvs.begin(); it != kvs.end(); it++) {
                node_element_count += it->bucket_element_count;
            }
            return node_element_count >= burst_threshold;
        }

        // find the key with keysize in its kvs
        T& access_kv_in_hashnode(const CharT* key, size_t keysize) {
            hash_node* targetNode = this;

            // TODO: burst if hash_node have too many element
            if (need_burst()) {
                std::vector<std::pair<std::string, T>> elements;
                for (auto it = kvs.begin(); it != kvs.end(); it++) {
                    // TODO: use the reference to reduce the extra copy
                    std::vector<std::pair<std::string, T>> bucket_element =
                        it->get_item_in_array_bucket();
                    elements.insert(elements.end(), bucket_element.begin(),
                                    bucket_element.end());
                }

                // update the targetNode
                targetNode = node_after_burst(elements, key, keysize);
            }

            size_t hashval = myTrie::hashRelative::hash<CharT>(key, keysize);
            std::cout << "bucket_id: " << hashval << std::endl;
            return targetNode->kvs[hashval].access_kv_in_bucket(key, keysize);
        }

        hash_node* node_after_burst(
            std::vector<std::pair<std::string, T>> kv_need_burst,
            const CharT* key, size_t keysize) {
            if (kv_need_burst.size() < burst_threshold) {
                hash_node* target = new hash_node(burst_threshold);
                        }
            std::map<CharT, std::vector<std::pair<std::string, T>>> newChilds;
            for (int i = 0; i != kv_need_burst.size(); i++) {
                std::pair<std::string, T> elem = kv_need_burst[i];
                newChilds[(elem.first)[0]].push_back(std::pair<std::string, T>(
                    elem.first.substr(1), elem.second));
            }
            for (auto it = newChilds.begin(); it != newChilds.end(); it++) {
                node_after_burst(it->second, key + 1, keysize - 1);
            }
            // TODO
        }

        // debug
        std::map<size_t, std::vector<std::pair<std::string, T>>>
        get_hash_nodes_buckets_info() {
            std::map<size_t, std::vector<std::pair<std::string, T>>> res;
            for (size_t i = 0; i != kvs.size(); i++) {
                res[i] = kvs[i].get_item_in_array_bucket();
            }
            return res;
        }
    };

    class array_bucket {
       public:
        CharT* arr_buffer;
        size_t bucket_element_count;
        size_t buffer_size;
        static const KeySizeT END_OF_BUCKET =
            std::numeric_limits<KeySizeT>::max();

        array_bucket() {
            bucket_element_count = 0;
            // inited array_bucket: |END_OF_BUCKET|
            arr_buffer = (CharT*)std::malloc(sizeof(END_OF_BUCKET));
            std::memcpy(arr_buffer, &END_OF_BUCKET, sizeof(END_OF_BUCKET));
            buffer_size = sizeof(END_OF_BUCKET);
        }

        bool is_end_of_bucket(const CharT* buffer) {
            return read_key_size(buffer) == END_OF_BUCKET;
        }

        size_t read_key_size(const CharT* buffer) {
            KeySizeT key_size;
            std::memcpy(&key_size, buffer, sizeof(KeySizeT));

            return (size_t)key_size;
        }

        // if found, return (entry_pos, true)
        // if not found, return (inserting_pos, false)
        std::pair<size_t, bool> find_in_bucket(const CharT* key,
                                               size_t keysize) {
            CharT* buffer_ptr = arr_buffer;
            size_t pos = 0;
            while (!is_end_of_bucket(buffer_ptr)) {
                size_t length = read_key_size(buffer_ptr);
                CharT* cmp_buffer_ptr = buffer_ptr + sizeof(KeySizeT);
                if (myTrie::hashRelative::keyEqual(key, keysize, cmp_buffer_ptr,
                                                   length)) {
                    return std::pair<size_t, bool>(pos, true);
                }
                // move ptr to next header, skip keysize, string, value
                buffer_ptr = buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                pos += sizeof(KeySizeT) + length + sizeof(T);
            }
            return std::pair<size_t, bool>(pos, false);
        }

        size_t cal_arrbuffer_newsize(size_t keysize) {
            // when add a new key we need extra |size|string|value|
            size_t new_buffer_size = buffer_size;
            new_buffer_size +=
                sizeof(KeySizeT) + keysize * sizeof(CharT) + sizeof(T);
            return new_buffer_size;
        }

        T& access_kv_in_bucket(const CharT* target, size_t keysize) {
            std::pair<size_t, bool> found = find_in_bucket(target, keysize);

            if (found.second) {
                // found the target, return the value reference
                return get_val_ref(found.first);

            } else {
                // not found, the buffer is full, need realloc
                size_t new_size = cal_arrbuffer_newsize(keysize);

                CharT* new_buffer =
                    (CharT*)(std::realloc(arr_buffer, new_size));
                if (new_buffer == nullptr) {
                    std::cout << "realloc failed!!!!!!!!!1" << std::endl;
                    exit(-1);
                }
                arr_buffer = new_buffer;
                buffer_size = new_size;

                T v = T{};
                append_impl(target, keysize, arr_buffer + found.first, v);
                bucket_element_count++;
                return get_val_ref(found.first);
            }
        }

        void append_impl(const CharT* target, size_t keysize,
                         CharT* buffer_append_pos, T& value) {
            // append key_size
            std::memcpy(buffer_append_pos, &keysize, sizeof(KeySizeT));
            buffer_append_pos += sizeof(KeySizeT);

            // append the string
            std::memcpy(buffer_append_pos, target, keysize * sizeof(CharT));
            buffer_append_pos += keysize;

            // append the value
            std::memcpy(buffer_append_pos, &value, sizeof(T));
            buffer_append_pos += sizeof(T);

            // append the whole buffer end
            const auto end_of_bucket = END_OF_BUCKET;
            std::memcpy(buffer_append_pos, &end_of_bucket,
                        sizeof(end_of_bucket));
        }

        T& get_val_ref(size_t val_pos) {
            CharT* entry_start_pos = arr_buffer + val_pos;
            size_t length = read_key_size(entry_start_pos);
            CharT* valAddress = entry_start_pos + sizeof(KeySizeT) + length;
            return *((T*)valAddress);
        }

        // debug
        void read_array_bucket() {
            std::cout << "\n\n-----------------------------\ncall "
                         "read_array_bucket\n";
            CharT* buf = arr_buffer;
            while (!is_end_of_bucket(buf)) {
                size_t length = read_key_size(buf);
                std::cout << *((KeySizeT*)buf) << "|";

                for (size_t i = 0; i != length; i++) {
                    std::cout << *(buf + sizeof(KeySizeT) + i) << "|";
                }
                // move ptr to next header, skip keysize, string, value
                buf = buf + sizeof(KeySizeT) + length;
                std::cout << (T)*buf << "|" << std::endl;
                buf = buf + +sizeof(T);
            }
            std::cout << "END read: |" << (int)*(buf) << "| finish "
                      << std::endl;
            std::cout << "-----------------------------\n\n";
        }

        std::vector<std::pair<std::string, T>> get_item_in_array_bucket() {
            // std::cout << "\n\n-----------------------------\ncall "
            //              "get_item_in_array_bucket\n";
            std::vector<std::pair<std::string, T>> res;
            CharT* buf = arr_buffer;
            while (!is_end_of_bucket(buf)) {
                size_t length = read_key_size(buf);

                char* temp = (char*)malloc(length + 1);
                std::memcpy(temp, buf + sizeof(KeySizeT), length);
                temp[length] = '\0';
                std::string item(temp);
                free(temp);

                // move ptr to next header, skip keysize, string, value
                buf = buf + sizeof(KeySizeT) + length;
                T v = *((T*)buf);

                res.push_back(std::pair<std::string, T>(item, v));
                buf = buf + sizeof(T);
                // std::cout << "get_item: " << item;
                // std::cout << " with size: " << item.size()
                //           << " with length: " << length
                //           << " with index4: " << (unsigned
                //           int)item[length]
                //           << "\n";
            }
            return res;
        }
    };

   public:
    anode* t_root;
    // class iterator {
    //    public:
    //     // trie_node* current_trie_node;
    //     hash_node* current_hash_node;

    //     size_t bucket_id;
    //     size_t buf_pos;
    // };
    htrie_map(size_t customized_burst_threshold = DEFAULT_BURST_THRESHOLD) {
        t_root = new hash_node(burst_threshold);
        burst_threshold = customized_burst_threshold;
    }

    T& operator[](std::string key) {
        return access_operator(key.data(), key.size());
    }

    T& operator[](const CharT* key) {
        return access_operator(key, std::strlen(key));
    }

    T& access_operator(const CharT* key, size_t key_size) {
        return find(key, key_size);
    }

    T& find(const CharT* key, size_t key_size) {
        anode* current_node = t_root;

        for (size_t pos = 0; pos < key_size; pos++) {
            // didn't reach the leaf
            // TODO: if a key end at a trie_node
            if (current_node->isTrieNode()) {
                current_node =
                    ((trie_node*)current_node)->findChildNode(key[pos]);
                // can't find
                if (current_node == nullptr) {
                }
            } else {
                // find a key in hash_node
                // if find: return the v reference
                // if notfind: write the key and return the v reference
                return ((hash_node*)current_node)
                    ->access_kv_in_hashnode(key + pos, key_size - pos);
            }
        }
    }
};

using std::cout;
using std::endl;
template <class CharT, class T>
void print_tree_struct(typename myTrie::htrie_map<CharT, T>::anode* root) {
    // print bucket
    if (root->isHashNode()) {
        std::map<size_t, std::vector<std::pair<std::string, T>>> buckets_info =
            ((typename myTrie::htrie_map<CharT, T>::hash_node*)root)
                ->get_hash_nodes_buckets_info();
        for (auto it = buckets_info.begin(); it != buckets_info.end(); it++) {
            cout << it->first << ":";
            std::vector<std::pair<std::string, T>> bucket_string = it->second;
            for (int i = 0; i != bucket_string.size(); i++) {
                cout << bucket_string[i].first << ",";
            }
            cout << "|" << endl;
        }
    } else if (root->isTrieNode()) {
        std::map<CharT, typename myTrie::htrie_map<CharT, T>::anode*> childs =
            ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                ->getChildsMap();
        for (auto it = childs.begin(); it != childs.end(); it++) {
            cout << it->first << "\t\t";
        }
        cout << "\nnext level:\n";
        for (auto it = childs.begin(); it != childs.end(); it++) {
            print_tree_struct<CharT, T>(it->second);
        }
    } else {
        cout << "node is not trie nor hash node\n";
        exit(0);
    }
}

template <class CharT, class T>
void print_htrie_map(htrie_map<CharT, T> hm,
                     std::vector<std::string> checklist) {
    // correctness check
    cout << "finding func check:\n";
    for (auto it = checklist.begin(); it != checklist.end(); it++) {
        cout << "search " << *it << " got " << hm[*it] << endl;
    }

    cout << "-------------------------\n";

    // node check
    print_tree_struct<CharT, T>(hm.t_root);
}
}  // namespace myTrie

int main() {
    using namespace std;
    using namespace myTrie;

    vector<string> tests;
    tests.push_back("abc");
    tests.push_back("abcd");
    tests.push_back("abcd");

    tests.push_back("abcde");

    htrie_map<char, int> hm;
    for (auto it = tests.begin(); it != tests.end(); it++)
        hm[*it] = it - tests.begin();
    std::cout << "------------------------\n";
    print_htrie_map<char, int>(hm, tests);
}