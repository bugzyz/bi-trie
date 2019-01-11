#include <stddef.h>
#include <climits>
#include <cstdint>

#include <cstdlib>
#include <vector>

#include <iostream>
#include <queue>
#include <sstream>
#include <stack>
#include "../util/hashFunc.hpp"

#include <fstream>
#include <ios>
#include <iostream>
#include <set>
#include <string>

#include <assert.h>

#define DEFAULT_Associativity 8
#define DEFAULT_Bucket_num 10

static uint32_t longest_string_size;

uint32_t recal_element_num_of_1st_char_counter = 0;
uint32_t burst_total_counter = 0;

// todo: wait to be deleted, just for recording the time that expand() cost
uint64_t expand_cost_time = 0;
uint64_t rehash_cost_time = 0;
uint64_t rehash_total_num = 0;

uint64_t shrink_total_time = 0;

uint64_t get_time() {
    struct timespec tp;
    /* POSIX.1-2008: Applications should use the clock_gettime() function
       instead of the obsolescent gettimeofday() function. */
    /* NOTE: The clock_gettime() function is only available on Linux.
       The mach_absolute_time() function is an alternative on OSX. */
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return ((tp.tv_sec * 1000 * 1000) + (tp.tv_nsec / 1000));
}

// configuration
// static size_t Associativity;
// static size_t Bucket_num;
// static size_t Max_bytes_per_kv;
// static double Burst_ratio;
// static size_t Max_slot_num;
// static size_t Max_loop;

size_t Associativity;
size_t Bucket_num;
size_t Max_bytes_per_kv;
size_t Max_slot_num;

using namespace std;

namespace myTrie {
// charT = char, T = value type, keysizeT = the type describe keysize
template <class CharT, class T, class KeySizeT = std::uint16_t>
class htrie_map {
   public:
    // DEFAULT_Associativity * DEFAULT_Bucket_num should be set to be greater
    // than 26 for 26 alaphbet and ", @ .etc.(the test in lubm40 shows that it
    // has 50 char species)
    enum class node_type : unsigned char {
        HASH_NODE,
        TRIE_NODE,
        MULTI_NODE,
    };

    class trie_node;
    class hash_node;
    class MultiTrieNode;
    class iterator;

    class anode {
       public:
        node_type _node_type;
        trie_node* parent;

        bool is_hash_node() { return _node_type == node_type::HASH_NODE; }
        bool is_trie_node() { return _node_type == node_type::TRIE_NODE; }
        bool is_multi_node() { return _node_type == node_type::MULTI_NODE; }

        void delete_me() {
            if (this->is_trie_node()) {
                std::map<CharT, trie_node*> childs =
                    ((trie_node*)this)->trie_node_childs;
                for (auto it = childs.begin(); it != childs.end(); it++) {
                    it->second->delete_me();
                    delete it->second;
                }
            } else {
                delete (hash_node*)this;
            }
        }
    };

    class multi_node : public anode {
       public:
        std::map<string, anode*> childs_;
        size_t string_keysize_;

        multi_node(size_t _string_keysize_)
            : string_keysize_(_string_keysize_) {
            anode::_node_type = node_type::MULTI_NODE;
            anode::parent = nullptr;
        }

        anode* find_child(const CharT* key, size_t key_size) {
            auto it = childs_.find(string(key, string_keysize_));
            if (it == childs_.end()) {
                return nullptr;
            } else {
                return it->second;
            }
        }
    };

    class trie_node : public anode {
       public:
        // store the suffix of hash_node or trie_node
        std::map<CharT, trie_node*> trie_node_childs;
        hash_node* hash_node_child;

        // prefix
        CharT* prefix;
        uint16_t prefix_len;

        // element end up here
        bool have_value;
        T value;

        trie_node(trie_node* p)
            : have_value(false),
              value(T()),
              hash_node_child(nullptr),
              prefix(nullptr),
              prefix_len(0) {
            anode::_node_type = node_type::TRIE_NODE;
            anode::parent = p;
        }

        iterator search_kv_in_trienode() {
            return iterator(have_value, value, this, 0, 0);
        }

        std::pair<bool, T> insert_kv_in_trienode(const CharT* key,
                                                 size_t key_size,
                                                 htrie_map<CharT, T>* hm, T v) {
            set_prefix(key, key_size);
            have_value = true;
            value = v;
            hm->set_v2k(v, this, 0, -1);
        }

        /*-------------------prefix relative-------------------*/
        size_t get_prefix(char* buf) {
            memcpy(buf, prefix, prefix_len);
            return prefix_len;
        }

        inline string get_prefix() {
            return prefix == nullptr ? string() : string(prefix, prefix_len);
        }

        void set_prefix(const CharT* key, size_t key_size) {
            char* cur_prefix = (char*)malloc(key_size);
            if (key_size != 0) {
                memcpy(cur_prefix, key, key_size);
                prefix = cur_prefix;
                prefix_len = key_size;
            } else {
                if (prefix != nullptr) free(prefix);
                prefix = nullptr;
                prefix_len = 0;
            }
        }

        void set_prefix(std::string& p) {
            uint16_t len = p.size();
            char* cur_prefix = (char*)malloc(p.size());
            if (len != 0) {
                memcpy(cur_prefix, p.data(), len);
                prefix = cur_prefix;
                prefix_len = len;
            } else {
                if (prefix != nullptr) free(prefix);
                prefix = nullptr;
                prefix_len = 0;
            }
        }

        ~trie_node() {
            std::map<CharT, trie_node*> empty;
            trie_node_childs.swap(empty);
            if (prefix != nullptr) free(prefix);
        }

        // finding target, if target doesn't exist, return nullptr
        trie_node* find_trie_node_child(CharT c) {
            auto found = trie_node_childs.find(c);
            if (found != trie_node_childs.end()) {
                trie_node* target = found->second;
                return target;
            } else {
                return nullptr;
            }
        }

        // finding target, if target doesn't exist, create new trie_node with
        // hash_node son and return new trie_node
        trie_node* find_trie_node_child(CharT c, const CharT* key,
                                        size_t key_size) {
            auto found = trie_node_childs.find(c);
            if (found != trie_node_childs.end()) {
                trie_node* target = found->second;
                return target;
            } else {
                trie_node* son_trie_node = new trie_node(this);
                this->add_trie_node_child(son_trie_node, c);
                son_trie_node->set_hash_node_child(
                    new hash_node(son_trie_node, string(key, key_size) + c));
                return son_trie_node;
            }
        }

        inline void set_hash_node_child(hash_node* node) {
            hash_node_child = node;
        }

        inline hash_node* get_hash_node_child() { return hash_node_child; }

        void add_trie_node_child(trie_node* node, CharT c) {
            trie_node_childs[c] = node;
            // clear the hash_node_child because a trie_node will only have
            // trie_node_childs map or have a single hash_node_child
            hash_node_child = nullptr;
        }

        pair<CharT, trie_node*> get_only_trie_node_child() {
            if (trie_node_childs.size() == 1) {
                auto it = trie_node_childs.begin();
                return pair<CharT, trie_node*>(it->first, it->second);
            }
            return pair<CharT, trie_node*>(CharT(), nullptr);
        }

        hash_node* get_only_hash_node_child() { return hash_node_child; }
    };

    class hash_node : public anode {
       public:
        class array_bucket;

        std::vector<array_bucket> kvs;

        size_t elem_num;
        size_t cur_page_id;

        int common_prefix_len;
        CharT* common_prefix;

        static const KeySizeT END_OF_BUCKET =
            std::numeric_limits<KeySizeT>::max();

        map<CharT, uint16_t> element_num_of_1st_char;

       public:
        inline static KeySizeT get_keysize(CharT* buffer) {
            return *(KeySizeT*)buffer;
        }

        inline static bool is_end_of_bucket(CharT* buffer) {
            return get_keysize(buffer) == END_OF_BUCKET;
        }

        inline static size_t get_need_size(size_t keysize) {
            return sizeof(KeySizeT) + keysize * sizeof(CharT) + sizeof(T);
        }

        inline CharT* get_first_key_pointer() { return common_prefix; }

        class array_bucket {
           public:
            CharT* arr_buffer;
            size_t buffer_size;

            array_bucket() {
                // inited array_bucket: |END_OF_BUCKET|
                arr_buffer = (CharT*)std::malloc(sizeof(KeySizeT));
                memcpy(arr_buffer, &END_OF_BUCKET, sizeof(END_OF_BUCKET));
                buffer_size = sizeof(END_OF_BUCKET);
            }

            // if found, return (entry_pos, true)
            // if not found, return (inserting_pos, false)
            std::pair<size_t, bool> find_in_bucket(const CharT* key,
                                                   size_t keysize) {
                CharT* buffer_ptr = arr_buffer;
                size_t pos = 0;
                while (!is_end_of_bucket(buffer_ptr)) {
                    size_t length = get_keysize(buffer_ptr);
                    CharT* cmp_buffer_ptr = buffer_ptr + sizeof(KeySizeT);
                    if (myTrie::hashRelative::keyEqual(cmp_buffer_ptr, length,
                                                       key, keysize)) {
                        return std::pair<size_t, bool>(pos, true);
                    }
                    // move ptr to next header, skip keysize, string, value
                    buffer_ptr =
                        buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                    pos += sizeof(KeySizeT) + length + sizeof(T);
                }
                return std::pair<size_t, bool>(pos, false);
            }

            void get_item_in_array_bucket(
                std::vector<std::pair<std::string, T>>& res) {
                CharT* buf = arr_buffer;
                while (!is_end_of_bucket(buf)) {
                    size_t length = get_keysize(buf);

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
                }
            }
        };

        void print_key_metas() {
            cout << "---------------\n";

            for (int i = 0; i != kvs.size(); i++) {
                vector<pair<string, T>> res;
                kvs[i].get_item_in_array_bucket(res);
                for (int ii = 0; ii != res.size(); ii++) {
                    auto it = res[ii];
                    cout << i << "." << ii << ": " << it.first << "="
                         << it.second << endl;
                }
                cout << "---\n";
            }
            cout << "---------------\n";
        }

       public:
        explicit hash_node(trie_node* p, string prefix)
            : elem_num(0),
              cur_page_id(0),
              common_prefix_len(INT_MAX),
              common_prefix(nullptr) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = p;
            if (p != nullptr) anode::parent->set_prefix(prefix);

            kvs = std::vector<array_bucket>(Bucket_num);
        }

        ~hash_node() {
            std::vector<array_bucket> empty;
            for (auto it = kvs.begin(); it != kvs.end(); it++) {
                free(it->arr_buffer);
            }
            kvs.swap(empty);
        }

        string get_prefix() {
            return anode::parent != nullptr ? anode::parent->get_prefix()
                                            : string();
        }

        inline char* get_tail_pointer(size_t bucketid, size_t pos) const {
            return kvs[bucketid].arr_buffer + pos;
        }

        inline T get_tail_v(char* buffer_ptr) {
            T v;
            std::memcpy(&v,
                        buffer_ptr + get_keysize(buffer_ptr) + sizeof(KeySizeT),
                        sizeof(T));
            return v;
        }

        void get_all_elements(std::map<std::string, T>& elements) {
            for (int i = 0; i != kvs.size(); i++) {
                CharT* buffer_ptr = kvs[i].arr_buffer;
                size_t pos = 0;
                while (!is_end_of_bucket(buffer_ptr)) {
                    size_t length = get_keysize(buffer_ptr);
                    CharT* cmp_buffer_ptr = buffer_ptr + sizeof(KeySizeT);

                    elements[string(buffer_ptr, length)] =
                        get_tail_v(buffer_ptr);

                    // move ptr to next header, skip keysize, string, value
                    buffer_ptr =
                        buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                    pos += sizeof(KeySizeT) + length + sizeof(T);
                }
            }
        }

        /*------------------ 3. bursting function------------------*/

        int cal_common_prefix_len(const char* s1, int cur_longest_prefix_len,
                                  const char* s2, int new_key_size) {
            if (cur_longest_prefix_len > new_key_size) {
                cur_longest_prefix_len = new_key_size;
            }
            for (int i = 0; i != cur_longest_prefix_len; i++) {
                if (s1[i] != s2[i]) {
                    return i;
                }
            }
            return cur_longest_prefix_len;
        }

        inline bool need_burst() const { return elem_num > Max_slot_num; }

        // To turn this(a hashnode) to n trie_node_childs of trie_node linking
        // their hashnode
        void burst(trie_node* p, htrie_map* hm, std::string prefix) {
            burst_total_counter++;
            // when the size of element_num_of_1st_char == 1, it means those
            // elements have a common prefix
            if (common_prefix_len != 0) {
                const char* first_key_p = get_first_key_pointer();

                // create the chain with several single trie_node
                // the number of node is common_prefix_len
                for (int i = 0; i != common_prefix_len; i++) {
                    if (p == nullptr) {
                        // bursting in a root hashnode
                        // the t_root is update to a empty trie_node
                        p = new trie_node(nullptr);
                        hm->setRoot(p);
                    }

                    trie_node* cur_trie_node = new trie_node(p);
                    p->add_trie_node_child(cur_trie_node, first_key_p[i]);
                    p = cur_trie_node;
                }

                // clear the element_num_of_1st_char
                element_num_of_1st_char.clear();

                // recalculate the capacity of hashnode we need
                for (int bucketid = 0; bucketid != kvs.size(); bucketid++) {
                    CharT* buffer_ptr = kvs[bucketid].arr_buffer;
                    size_t pos = 0;
                    while (!is_end_of_bucket(buffer_ptr)) {
                        size_t length = get_keysize(buffer_ptr);
                        char* key = buffer_ptr + sizeof(KeySizeT);

                        element_num_of_1st_char[key[common_prefix_len]]++;

                        // move ptr to next header, skip keysize, string, value
                        buffer_ptr =
                            buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                        pos += sizeof(KeySizeT) + length + sizeof(T);
                    }
                }

                // debug
                recal_element_num_of_1st_char_counter++;

                // update prefix to (prior prefix + common chain prefix)
                prefix = prefix + string(first_key_p, common_prefix_len);
            }

            map<char, hash_node*> hnode_set;
            // create several hashnode based on the number of the elements
            // that start with the same char
            for (auto it = element_num_of_1st_char.begin();
                 it != element_num_of_1st_char.end(); it++) {
                if (p == nullptr) {
                    // bursting in a root hashnode
                    // the t_root is update to a empty trie_node
                    p = new trie_node(nullptr);
                    hm->setRoot(p);
                }
                trie_node* cur_trie_node = new trie_node(p);
                p->add_trie_node_child(cur_trie_node, it->first);

                hash_node* hnode =
                    new hash_node(cur_trie_node, prefix + it->first);
                cur_trie_node->set_hash_node_child(hnode);

                hnode_set[it->first] = hnode;
            }

            // a map for the hashnode that insert fail
            // map<char, pair<hash_node*, map<string, T>>> burst_again_list;

            // insert the elements with same first char after common_prefix_len
            for (int bucketid = 0; bucketid != kvs.size(); bucketid++) {
                CharT* buffer_ptr = kvs[bucketid].arr_buffer;
                size_t pos = 0;
                while (!is_end_of_bucket(buffer_ptr)) {
                    size_t length = get_keysize(buffer_ptr);
                    char* key = buffer_ptr + sizeof(KeySizeT);

                    size_t length_left = length - common_prefix_len - 1;

                    T v = get_tail_v(buffer_ptr);
                    hash_node* hnode = hnode_set[key[common_prefix_len]];

                    // rarely, insert in trie_node
                    if (length_left == 0) {
                        string new_prefix = prefix + key[common_prefix_len];

                        hnode->anode::parent->insert_kv_in_trienode(
                            new_prefix.data(), new_prefix.size(), hm, v);

                        // move ptr to next header, skip keysize, string, value
                        buffer_ptr =
                            buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                        pos += sizeof(KeySizeT) + length + sizeof(T);
                        continue;
                    }

                    // normally, insert in hash_node
                    iterator target_it = hnode->search_kv_in_hashnode(
                        key + common_prefix_len + 1, length_left);
                    std::pair<bool, T> res = target_it.insert_hashnode(
                        key + common_prefix_len + 1, length_left, hm, v);

                    // move ptr to next header, skip keysize, string, value
                    buffer_ptr =
                        buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                    pos += sizeof(KeySizeT) + length + sizeof(T);
                }
            }
            return;
        }

        /*----------------searching in hash_node----------------*/

        // return <found?, iterator>
        std::pair<bool, iterator> find_in_bucket(size_t bucketid,
                                                 const CharT* key,
                                                 size_t keysize) {
            CharT* buffer_ptr = kvs[bucketid].arr_buffer;
            size_t pos = 0;
            while (!is_end_of_bucket(buffer_ptr)) {
                size_t length = get_keysize(buffer_ptr);
                CharT* cmp_buffer_ptr = buffer_ptr + sizeof(KeySizeT);
                if (myTrie::hashRelative::keyEqual(cmp_buffer_ptr, length, key,
                                                   keysize)) {
                    return std::pair<bool, iterator>(
                        true, iterator(true, get_tail_v(buffer_ptr), this,
                                       bucketid, pos));
                }
                // move ptr to next header, skip keysize, string, value
                buffer_ptr = buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
                pos += sizeof(KeySizeT) + length + sizeof(T);
            }
            // don't find in this bucket, so we return the iterator with current
            // bucketid and pos if user want to add new element
            return std::pair<bool, iterator>(
                false, iterator(false, T(), this, bucketid, pos));
        }

        iterator search_kv_in_hashnode(const CharT* key, size_t keysize) {
            // if found the existed target in bucket1 or bucket2, just
            // return the iterator for being modified or read
            size_t bucketId1 =
                myTrie::hashRelative::hash(key, keysize, 1) % Bucket_num;
            std::pair<bool, iterator> res1 =
                find_in_bucket(bucketId1, key, keysize);

            return res1.second;
        }

        /*----------------inserting in hash_node----------------*/

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
            std::memcpy(buffer_append_pos, &END_OF_BUCKET,
                        sizeof(END_OF_BUCKET));
        }

        // return: inserting pointer
        char* alloc_insert_space(size_t bucketid, size_t keysize) {
            // calculate the new size we gonna realloc
            size_t new_size =
                kvs[bucketid].buffer_size + get_need_size(keysize);
            size_t ret_pos = kvs[bucketid].buffer_size - sizeof(KeySizeT);

            // update arr_buffer and buffer_size
            kvs[bucketid].arr_buffer =
                (char*)realloc(kvs[bucketid].arr_buffer, new_size);
            kvs[bucketid].buffer_size = new_size;

            return kvs[bucketid].arr_buffer + ret_pos;
        }

        std::pair<bool, T> insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, size_t bucketid,
                                                 uint16_t pos,
                                                 std::string prefix) {
            if (elem_num == 0) {
                common_prefix = (CharT*)malloc(keysize * sizeof(CharT));
                memcpy(common_prefix, key, keysize);
            }

            // allocate new page or alloc more space in old page
            char* inserting_pointer = alloc_insert_space(bucketid, keysize);

            append_impl(key, keysize, inserting_pointer, v);

            // set v2k
            hm->set_v2k(v, this, bucketid, pos);
            elem_num++;

            // update the element_num_of_1st_char
            element_num_of_1st_char[*key]++;

            // update the common_prefix_len
            int cur_com_prefix_len = cal_common_prefix_len(
                get_first_key_pointer(), common_prefix_len, key, keysize);
            if (common_prefix_len > cur_com_prefix_len) {
                common_prefix_len = cur_com_prefix_len;
            }

            if (need_burst()) {
                burst(this->anode::parent, hm, prefix);

                delete this;
            }

            return std::pair<bool, T>(true, v);
        }
    };

    class SearchPoint {
       public:
       public:
        anode* node;
        int32_t bucketid;
        int32_t pos;

        SearchPoint() : node(nullptr), bucketid(-1), pos(-1) {}
        SearchPoint(anode* n, int32_t bi, int32_t p)
            : node(n), bucketid(bi), pos(p) {}

        std::string get_string() {
            if (node == nullptr) return string();
            // if the node is trie_node, just return the prefix on node
            if (node->is_trie_node()) {
                return ((trie_node*)node)->get_prefix();
            }

            // get the parent char chain
            trie_node* cur_node = ((hash_node*)node)->anode::parent;
            char* buf = (char*)malloc(longest_string_size);

            size_t len = cur_node != nullptr ? cur_node->get_prefix(buf) : 0;

            // get tail
            if (pos != -1) {
                char* buf_ptr =
                    ((hash_node*)node)->get_tail_pointer(bucketid, pos);
                size_t keysize = hash_node::get_keysize(buf_ptr);

                memcpy(buf + len, buf_ptr + sizeof(KeySizeT), keysize);
                len += keysize;
            }
            string res = string(buf, len);
            free(buf);
            return res;
        }
    };

   public:
    anode* t_root;
    std::map<T, SearchPoint> v2k;
    htrie_map(size_t customized_associativity = DEFAULT_Associativity,
              size_t customized_bucket_count = DEFAULT_Bucket_num)
        : t_root(nullptr) {
        std::cout << "SET UP TESSIL-HAT-TRIE MAP\n";
        cout << "NO GROW\n";

        Associativity = customized_associativity;
        Bucket_num = customized_bucket_count;

        Max_slot_num = Associativity * Bucket_num;

        t_root = new hash_node(nullptr, string());
    }

    void set_searchPoint_index(T v, int index) { v2k[v].set_index(index); }

    void set_v2k(T v, anode* node, int32_t bucketid, int32_t pos) {
        v2k[v] = SearchPoint(node, bucketid, pos);
    }

    void setRoot(anode* node) { t_root = node; }

    class iterator {
       public:
        bool found;
        const T v;
        anode* target_node;
        size_t bucketid;
        uint16_t pos;

        iterator(bool f, T vv, anode* hnode, size_t bid, int p)
            : found(f), v(vv), target_node(hnode), bucketid(bid), pos(p) {}

        std::pair<bool, T> insert_hashnode(const CharT* key, size_t key_size,
                                           htrie_map<CharT, T>* hm, T v) {
            return ((hash_node*)target_node)
                ->insert_kv_in_hashnode(
                    key, key_size, hm, v, bucketid, pos,
                    ((hash_node*)target_node)->get_prefix());
        }

        std::pair<bool, T> insert_trienode(const CharT* key, size_t key_size,
                                           htrie_map<CharT, T>* hm, T v) {
            return ((trie_node*)target_node)
                ->insert_kv_in_trienode(key, key_size, hm, v);
        }
    };

    std::pair<bool, T> access_kv_in_htrie_map(const CharT* key, size_t key_size,
                                              T v, bool findMode) {
        // update longest_string_size
        longest_string_size =
            longest_string_size > key_size ? longest_string_size : key_size;

        anode* current_node = t_root;

        for (size_t pos = 0; pos < key_size; pos++) {
            switch (current_node->anode::_node_type) {
                case node_type::MULTI_NODE: {
                    if (!findMode) {
                    } else {
                        size_t jump_pos =
                            ((multi_node*)current_node)->string_keysize_;
                        current_node =
                            ((multi_node*)current_node)
                                ->find_child(key + pos, key_size - pos);
                        if (current_node == nullptr) {
                            return std::pair<bool, T>(false, T());
                        }
                        pos += jump_pos - 1;

                        if (current_node->is_trie_node() &&
                            ((trie_node*)current_node)->get_hash_node_child() !=
                                nullptr &&
                            pos != key_size - 1) {
                            current_node = ((trie_node*)current_node)
                                               ->get_hash_node_child();
                        }
                    }
                } break;
                case node_type::TRIE_NODE: {
                    if (!findMode) {
                        // return the hitted trie_node* or create a new
                        // trie_node with a hash_node son
                        current_node =
                            ((trie_node*)current_node)
                                ->find_trie_node_child(key[pos], key, pos);
                    } else {
                        // return the hitted trie_node* or nullptr if not
                        // found
                        current_node = ((trie_node*)current_node)
                                           ->find_trie_node_child(key[pos]);
                        // only in the findMode==true can cause the
                        // current_node to be nullptr
                        if (current_node == nullptr) {
                            return std::pair<bool, T>(false, T());
                        }
                    }

                    if (((trie_node*)current_node)->get_hash_node_child() !=
                            nullptr &&
                        pos != key_size - 1) {
                        current_node =
                            ((trie_node*)current_node)->get_hash_node_child();
                    }
                } break;
                case node_type::HASH_NODE: {
                    iterator it =
                        ((hash_node*)current_node)
                            ->search_kv_in_hashnode(key + pos, key_size - pos);
                    if (findMode) {
                        return std::pair<bool, T>(it.found, it.v);
                    } else {
                        pair<bool, T> res = it.insert_hashnode(
                            key + pos, key_size - pos, this, v);
                        if (res.first == false) {
                            // if the insert failed, we burst the
                            // target_hashnode and retry insertion
                            hash_node* hnode_burst_needed =
                                (hash_node*)current_node;

                            hnode_burst_needed->burst(
                                hnode_burst_needed->anode::parent, this,
                                string(key, pos));

                            delete hnode_burst_needed;
                            return access_kv_in_htrie_map(key, key_size, v,
                                                          false);
                        }
                        return res;
                    }
                } break;
                default:
                    cout << "wrong type!";
                    exit(0);
            }
        }

        // find a key in trie_node's only value
        iterator it = ((trie_node*)current_node)->search_kv_in_trienode();

        if (findMode) {
            return std::pair<bool, T>(it.found, it.v);
        } else {
            return it.insert_trienode(key, key_size, this, v);
        }
    }

    /*---------------external accessing interface-------------------*/

    // access element
    T searchByKey(std::string key) {
        return access_kv_in_htrie_map(key.data(), key.size(), T(), true).second;
    }

    std::string searchByValue(T v) { return v2k[v].get_string(); }

    // find operation
    std::pair<bool, T> findByKey(std::string key) {
        return access_kv_in_htrie_map(key.data(), key.size(), T(), true);
    }

    std::pair<bool, std::string> findByValue(T v) {
        if (v2k.find(v) == v2k.end()) {
            return std::pair<bool, T>(false, string());
        } else {
            return std::pair<bool, T>(true, v2k[v].get_string());
        }
    }

    std::pair<bool, T> insertKV(std::string key, T v) {
        return access_kv_in_htrie_map(key.data(), key.size(), v, false);
    }

    /*---------------external cleaning interface-------------------*/

    anode* shrink_node(anode* node) {
        if (node->is_trie_node()) {
            trie_node* cur_node = (trie_node*)node;

            if (cur_node->get_only_hash_node_child() != nullptr) {
                return cur_node;
            }

            vector<pair<string, anode*>> traverse_save(
                cur_node->trie_node_childs.size());
            vector<pair<CharT, trie_node*>> next_layer;

            for (auto it = cur_node->trie_node_childs.begin();
                 it != cur_node->trie_node_childs.end(); it++) {
                next_layer.push_back(
                    pair<CharT, trie_node*>(it->first, it->second));
            }

            bool allow_next_layer = true;
            size_t string_keysize = 0;
            do {
                for (int i = 0; i != next_layer.size(); i++) {
                    CharT c = next_layer[i].first;
                    trie_node* next_layer_trie_node = next_layer[i].second;

                    // current key_string add the next_layer's char
                    string new_key_string = traverse_save[i].first + c;

                    traverse_save[i].first = new_key_string;
                    traverse_save[i].second = next_layer_trie_node;

                    pair<CharT, trie_node*> next_pair;
                    if (next_layer_trie_node->get_only_hash_node_child() ==
                            nullptr &&
                        next_layer_trie_node->get_only_trie_node_child()
                                .second == nullptr) {
                        // have several children
                        allow_next_layer = false;
                    } else {
                        // only have a hash_node child or only have a
                        // trie_node child
                        if (next_layer_trie_node->get_only_hash_node_child() !=
                            nullptr) {
                            next_layer[i] =
                                pair<CharT, trie_node*>(CharT(), nullptr);
                            // stop at this layer
                            allow_next_layer = false;
                            traverse_save[i].second = next_layer_trie_node;
                        } else {
                            next_pair = next_layer_trie_node
                                            ->get_only_trie_node_child();
                            next_layer[i] = next_pair;
                        }
                    }
                }
                if (allow_next_layer) {
                    for (int i = 0; i != next_layer.size(); i++) {
                        delete (next_layer[i].second)->anode::parent;
                    }
                }
                string_keysize++;
            } while (allow_next_layer);

            // construct the target multi_node
            multi_node* target_node = new multi_node(string_keysize);

            for (int i = 0; i != traverse_save.size(); i++) {
                anode* res = shrink_node(traverse_save[i].second);
                target_node->childs_[traverse_save[i].first] = res;
            }
            return target_node;
        } else if (node->is_hash_node()) {
            return node;
        } else {
            cout << "program are in a unexpected branch\n";
            assert(false);
            exit(0);
        }
    }

    void shrink() {
        cout << "Shrinking\n";
        uint64_t sta = get_time();
        t_root = shrink_node(t_root);
        uint64_t end = get_time();
        shrink_total_time = end - sta;
    }

    void deleteMyself() {
        map<T, SearchPoint> empty;
        v2k.swap(empty);
        t_root->delete_me();
    }
};  // namespace myTrie

}  // namespace myTrie