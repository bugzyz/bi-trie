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
#define DEFAULT_Max_bytes_per_kv 1000
#define DEFAULT_Burst_ratio 0.75

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
double Burst_ratio;
size_t Max_slot_num;
size_t Max_loop;

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
            hm->set_v2k(v, this, -1);
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
        class slot;

        slot* key_metas;
        vector<std::pair<char*, size_t>> pages;
        size_t elem_num;
        size_t cur_page_id;

        int common_prefix_len;

        map<CharT, uint16_t> element_num_of_1st_char;

       public:
        // debug function
        void print_slot(int i, int j) {
            slot* s = key_metas + i * Associativity + j;
            string str = string(get_tail_pointer(s), s->length);
            T v = get_tail_v(s);
            cout << i * Associativity + j << ":" << s->length << "," << s->pos
                 << "," << s->page_id << "," << str << "=" << v << "\n";
        }

        void print_key_metas() {
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != Associativity; j++) {
                    print_slot(i, j);
                }
                cout << "---\n";
            }
        }

       public:
        struct slot {
            KeySizeT length;
            uint16_t pos;
            uint16_t page_id;

            bool isEmpty() { return length == 0; }

            slot(KeySizeT l, size_t p, size_t pi)
                : length(l), pos(p), page_id(pi) {}

            void set_slot(KeySizeT l, size_t p, size_t pi) {
                length = l;
                pos = p;
                page_id = pi;
            }

            void set_slot(slot* s) {
                length = s->length;
                pos = s->pos;
                page_id = s->page_id;
            }
        };

       public:
        explicit hash_node(trie_node* p, string prefix)
            : elem_num(0), cur_page_id(0), common_prefix_len(INT_MAX) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = p;
            if (p != nullptr) anode::parent->set_prefix(prefix);

            key_metas =
                (slot*)malloc(Associativity * Bucket_num * sizeof(slot));

            // init key space
            for (int i = 0; i != Associativity * Bucket_num; i++) {
                key_metas[i].length = 0;
                key_metas[i].pos = 0;
                key_metas[i].page_id = 0;
            }

            char* page = (char*)malloc(Max_bytes_per_kv);
            pages.push_back(std::pair<char*, size_t>(page, 0));
        }

        ~hash_node() {
            free(key_metas);
            for (int i = 0; i != pages.size(); i++) free(pages[i].first);
            vector<std::pair<char*, size_t>> empty;
            pages.swap(empty);
        }

        string get_prefix() {
            return anode::parent != nullptr ? anode::parent->get_prefix()
                                            : string();
        }

        inline char* get_tail_pointer(slot* s) const {
            return pages[s->page_id].first + s->pos;
        }

        inline slot* get_slot(size_t bucketid, size_t slotid) {
            return key_metas + bucketid * Associativity + slotid;
        }

        inline slot* get_slot(int index) { return key_metas + index; }

        inline int get_index(slot* s) { return s - key_metas; }

        void get_tail_str_v(std::map<std::string, T>& elements, slot* s) {
            char* tail_pointer = get_tail_pointer(s);
            std::string res(tail_pointer, s->length);
            std::memcpy(&elements[res], tail_pointer + s->length, sizeof(T));
        }

        T get_tail_v(slot* s) {
            T v;
            std::memcpy(&v, get_tail_pointer(s) + s->length, sizeof(T));
            return v;
        }

        void get_all_elements(std::map<std::string, T>& elements) {
            for (size_t i = 0; i != Bucket_num; i++) {
                for (size_t j = 0; j != Associativity; j++) {
                    slot& cur_slot = key_metas[i * Associativity + j];
                    if (cur_slot.isEmpty()) break;
                    get_tail_str_v(elements, &cur_slot);
                }
            }
        }

        /*---------function that changed key_metas layout---------*/

        /*------------------ 0. helper function------------------*/

        void apply_the_changed_searchPoint(map<T, int>& searchPoints,
                                           htrie_map<CharT, T>* hm) {
            for (auto it = searchPoints.begin(); it != searchPoints.end(); it++)
                hm->set_searchPoint_index(it->first, it->second);
        }

        /*------------------ 2. rehashing function------------------*/

        inline slot* previous_dst_slot_in_same_bucket(slot* s) {
            size_t slotid = (s - key_metas) % Associativity;
            if (slotid == 0) {
                return nullptr;
            } else
                return s - 1;
        }

        // return another possible bucketid that the slot *s can be
        inline size_t get_another_bucketid(slot* s, size_t current_bucketid) {
            size_t bucketid1 =
                myTrie::hashRelative::hash(get_tail_pointer(s), s->length, 1) %
                Bucket_num;
            size_t bucketid2 =
                myTrie::hashRelative::hash(get_tail_pointer(s), s->length, 2) %
                Bucket_num;
            return current_bucketid == bucketid1 ? bucketid2 : bucketid1;
        }

        int cuckoo_hash(size_t bucketid, htrie_map<CharT, T>* hm) {
            rehash_total_num++;
            uint64_t sta = get_time();
            // bucket_list records the mapping of bucket_id=last_empty_slot_id
            std::map<size_t, size_t> bucket_list;
            for (size_t bn = 0; bn != Bucket_num; bn++) {
                bucket_list[bn] = Associativity;
                for (int sn = 0; sn != Associativity; sn++) {
                    if (key_metas[bn * Associativity + sn].isEmpty()) {
                        bucket_list[bn] = sn;
                        break;
                    }
                }
            }
            // current bucket is definitely full
            // just pick the last slot to kick
            int ret_slot_id = Associativity - 1;

            size_t kicked_slot_id = -1;
            for (int i = 0; i != Associativity; i++) {
                slot* s = get_slot(bucketid, i);
                size_t bkid = get_another_bucketid(s, bucketid);
                if (bkid != bucketid) {
                    kicked_slot_id = i;
                }
            }
            if (kicked_slot_id == -1) {
                uint64_t end = get_time();
                rehash_cost_time += end - sta;

                return -1;
            }
            // set up the backup for recovery if the rehash fails
            char* key_metas_backup =
                (char*)malloc(Bucket_num * Associativity * sizeof(slot));
            memcpy(key_metas_backup, key_metas,
                   Bucket_num * Associativity * sizeof(slot));

            ret_slot_id = kicked_slot_id;

            size_t current_bucket_id = bucketid;

            slot src_slot = slot(0, 0, 0);
            slot* dst_slot = get_slot(current_bucket_id, kicked_slot_id);

            size_t rehash_count = 0;

            size_t last_current_bucketid = 0;
            size_t last_bucketid_kick_to = 0;

            map<T, int> searchPoint_wait_2_be_update;
            do {
                /*
                    src(a,b,c)
                    cur_bucket: |x      |x      |x      |dst(d,e,f)|
                    kk2_bucket: |x      |x      |x      |x         |
                */
                // calculate the destination
                size_t bucketid_kick_to =
                    get_another_bucketid(dst_slot, current_bucket_id);

                // if the slot can only place in one bucket, we change the
                // dst_slot
                // if the cuckoo hash kick as a circle, we change the dst_slot
                // to try to break the circle
                if (bucketid_kick_to == current_bucket_id ||
                    (last_bucketid_kick_to == current_bucket_id &&
                     last_current_bucketid == bucketid_kick_to)) {
                    dst_slot = previous_dst_slot_in_same_bucket(dst_slot);
                    rehash_count++;
                    if (dst_slot == nullptr) {
                        // recover the key_metas
                        memcpy(key_metas, key_metas_backup,
                               Bucket_num * Associativity * sizeof(slot));
                        free(key_metas_backup);
                        uint64_t end = get_time();
                        rehash_cost_time += end - sta;
                        return -1;
                    }
                    continue;
                }

                /*
                    src(a,b,c)
                    temp: d,e,f
                    cur_bucket: |x      |x      |x      |dst(d,e,f)|
                    kk2_bucket: |x      |x      |x      |x         |
                */
                KeySizeT temp_length = dst_slot->length;
                size_t temp_pos = dst_slot->pos;
                size_t temp_page_id = dst_slot->page_id;

                // if dst_slot is empty, it means now the dst_slot is the first
                // place we clear for the target slot
                if (dst_slot->isEmpty()) {
                    // recover the key_metas
                    memcpy(key_metas, key_metas_backup,
                           Bucket_num * Associativity * sizeof(slot));
                    free(key_metas_backup);
                    uint64_t end = get_time();
                    rehash_cost_time += end - sta;
                    return -1;
                }

                /*
                    src(a,b,c)
                    temp: d,e,f
                    cur_bucket: |x      |x      |x      |dst(a,b,c)|
                    kk2_bucket: |x      |x      |x      |x         |
                */
                dst_slot->set_slot(src_slot.length, src_slot.pos,
                                   src_slot.page_id);
                if (!dst_slot->isEmpty())
                    searchPoint_wait_2_be_update[get_tail_v(dst_slot)] =
                        get_index(dst_slot);

                // if the destination bucket isn't full, just fill the empty
                // slot and return
                if (bucket_list[bucketid_kick_to] != Associativity) {
                    /* kick2bucket is have empty slot
                        src(a,b,c)
                        temp: d,e,f
                        cur_bucket: |x      |x      |x      |dst(a,b,c)|
                        kk2_bucket: |x      |x      |x      |0         |
                    */
                    // update dst_slot to the dest bucket's first empty slot
                    /*
                        src(a,b,c)
                        temp: d,e,f
                        cur_bucket: |x      |x      |x      |x(a,b,c)  |
                        kk2_bucket: |x      |x      |x      |0-dst     |
                    */
                    dst_slot = get_slot(bucketid_kick_to,
                                        bucket_list[bucketid_kick_to]);

                    /*
                        src(a,b,c)
                        temp: d,e,f
                        cur_bucket: |x      |x      |x      |x(a,b,c)  |
                        kk2_bucket: |x      |x      |x      |dst(d,e,f)|
                    */
                    dst_slot->set_slot(temp_length, temp_pos, temp_page_id);
                    searchPoint_wait_2_be_update[get_tail_v(dst_slot)] =
                        get_index(dst_slot);
                    apply_the_changed_searchPoint(searchPoint_wait_2_be_update,
                                                  hm);
                    free(key_metas_backup);
                    uint64_t end = get_time();
                    rehash_cost_time += end - sta;
                    return ret_slot_id;
                }

                /* kick2bucket is full
                    src(a,b,c)
                    temp: d,e,f
                    cur_bucket: |x      |x      |x      |dst(a,b,c)|
                    kk2_bucket: |x      |x      |x      |x         |
                */

                /* kick2bucket is full
                    src(a,b,c)
                    temp: d,e,f
                    cur_bucket: |x      |x      |x      |x(a,b,c)  |
                    kk2_bucket: |x      |x      |x      |x-dst     |
                */
                // update dst_slot to the dest bucket's last empty slot
                dst_slot = get_slot(bucketid_kick_to, Associativity - 1);
                /*
                     src(d,e,f)
                     cur_bucket: |x      |x      |x      |x(a,b,c)  |
                     kk2_bucket: |x      |x      |x      |x-dst     |
                */
                src_slot.set_slot(temp_length, temp_pos, temp_page_id);
                /*
                     src(d,e,f)
                     (kk2_bucket)
                     cur_bucket: |x      |x      |x      |x-dst     |
                */
                last_bucketid_kick_to = bucketid_kick_to;
                last_current_bucketid = current_bucket_id;
                current_bucket_id = bucketid_kick_to;

                rehash_count++;
            } while (rehash_count != Max_loop);
            // recover the key_metas
            memcpy(key_metas, key_metas_backup,
                   Bucket_num * Associativity * sizeof(slot));
            free(key_metas_backup);
            uint64_t end = get_time();
            rehash_cost_time += end - sta;
            return -1;
        }

        /*------------------ 3. bursting function------------------*/

        inline const char* get_first_key_pointer() { return pages[0].first; }

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

        inline size_t round_up_2_next_power_2(size_t x) {
            x--;
            x |= x >> 1;   // handle  2 bit numbers
            x |= x >> 2;   // handle  4 bit numbers
            x |= x >> 4;   // handle  8 bit numbers
            x |= x >> 8;   // handle 16 bit numbers
            x |= x >> 16;  // handle 32 bit numbers
            x++;

            return x;
        }

        inline bool need_burst() const {
            return elem_num >= Max_slot_num * Burst_ratio;
        }

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
                for (int i = 0; i != Bucket_num; i++) {
                    for (int j = 0; j != Associativity; j++) {
                        slot* s = get_slot(i, j);
                        if (s->isEmpty()) break;

                        char* key = get_tail_pointer(s);
                        element_num_of_1st_char[key[common_prefix_len]]++;
                    }
                }

                // debug
                recal_element_num_of_1st_char_counter++;

                // update prefix to (prior prefix + common chain prefix)
                prefix = prefix + string(first_key_p, common_prefix_len);
            }

            map<char, hash_node*> hnode_set;
            // create several hashnode based on the number of the elements that
            // start with the same char
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
            map<char, pair<hash_node*, map<string, T>>> burst_again_list;

            // insert the elements with same first char after common_prefix_len
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != Associativity; j++) {
                    slot* s = get_slot(i, j);
                    if (s->isEmpty()) break;

                    char* key = get_tail_pointer(s);

                    size_t length_left = s->length - common_prefix_len - 1;

                    T v = get_tail_v(s);
                    hash_node* hnode = hnode_set[key[common_prefix_len]];

                    if (hnode == nullptr) {
                        map<string, T>& temp_elements =
                            burst_again_list[key[common_prefix_len]].second;
                        temp_elements[string(key + common_prefix_len + 1,
                                             length_left)] = v;
                        continue;
                    }

                    // rarely, insert in trie_node
                    if (length_left == 0) {
                        string new_prefix = prefix + key[common_prefix_len];

                        hnode->anode::parent->insert_kv_in_trienode(
                            new_prefix.data(), new_prefix.size(), hm, v);
                        continue;
                    }

                    // normally, insert in hash_node
                    iterator target_it = hnode->search_kv_in_hashnode(
                        key + common_prefix_len + 1, length_left);
                    std::pair<bool, T> res = target_it.insert_hashnode(
                        key + common_prefix_len + 1, length_left, hm, v);

                    // if insert failed, it need burst again
                    if (res.first == false) {
                        hnode_set[key[common_prefix_len]] = nullptr;

                        // get the element already inserted
                        map<string, T> temp_elements;
                        hnode->get_all_elements(temp_elements);

                        // and current key=value
                        temp_elements[string(key + common_prefix_len + 1,
                                             length_left)] = v;

                        burst_again_list[key[common_prefix_len]] =
                            pair<hash_node*, map<string, T>>(hnode,
                                                             temp_elements);
                    }
                }
            }

            // check if there is some failure when insert in hash_node
            for (auto it = burst_again_list.begin();
                 it != burst_again_list.end(); it++) {
                // If burst() fail at this prefix char, we use the old way to
                // burst
                string burst_again_prefix = prefix + it->first;
                pair<hash_node*, map<string, T>>& temp_pair = it->second;
                hash_node* burst_hnode = temp_pair.first;
                burst_hnode->burst_by_elements(temp_pair.second,
                                               burst_hnode->anode::parent, hm,
                                               burst_again_prefix);
            }
            return;
        }

        // optional burst function, time-consuming, but work fine
        // To turn this(a hashnode) to n trie_node_childs of trie_node linking
        // their hashnode
        void burst_by_elements(std::map<std::string, T>& elements, trie_node* p,
                               htrie_map* hm, std::string prefix) {
            // calculate the prefix len first
            auto it = elements.begin();
            string first_string = it->first;
            const char* first_key_p = first_string.data();
            int common_prefix_len = first_string.size();

            // those hash_node that without inserting all elements will have a
            // wrong common_prefix_len, so here we recalculate the
            // common_prefix_len
            for (it++; it != elements.end(); it++) {
                const string& cur_string = it->first;
                int cur_common_prefix =
                    cal_common_prefix_len(first_key_p, common_prefix_len,
                                          cur_string.data(), cur_string.size());
                if (common_prefix_len > cur_common_prefix) {
                    common_prefix_len = cur_common_prefix;
                }
            }

            // create the chain with several single trie_node
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

            // after now, the subsequent key are different at first
            // char
            std::map<CharT, std::map<std::string, T>> preprocElements;
            for (auto it = elements.begin(); it != elements.end(); it++) {
                const string& cur_string = it->first;
                preprocElements[cur_string[common_prefix_len]]
                               [cur_string.substr(common_prefix_len + 1)] =
                                   it->second;
            }

            // update prefix to (prior prefix + common chain prefix)
            prefix = prefix + string(first_key_p, common_prefix_len);

            // deal with the element with several different head char
            for (auto it = preprocElements.begin(); it != preprocElements.end();
                 it++) {
                if (p == nullptr) {
                    // bursting in a root hashnode
                    // the t_root is update to a empty trie_node
                    p = new trie_node(nullptr);
                    hm->setRoot(p);
                }
                trie_node* cur_trie_node = new trie_node(p);
                p->add_trie_node_child(cur_trie_node, it->first);

                std::map<std::string, T>& curKV = it->second;

                hash_node* hnode =
                    new hash_node(cur_trie_node, prefix + it->first);
                cur_trie_node->set_hash_node_child(hnode);

                bool stop_insert_and_burst = false;
                for (auto itt = curKV.begin(); itt != curKV.end(); itt++) {
                    const string& temp = itt->first;

                    if (temp.size() == 0) {
                        cur_trie_node->insert_kv_in_trienode(
                            (prefix + it->first).data(), prefix.size() + 1, hm,
                            itt->second);
                        continue;
                    }

                    iterator target_it =
                        hnode->search_kv_in_hashnode(temp.data(), temp.size());
                    std::pair<bool, T> res = target_it.insert_hashnode(
                        temp.data(), temp.size(), hm, itt->second);
                    // if insert failed, it need burst
                    if (res.first == false) {
                        stop_insert_and_burst = true;
                        break;
                    }
                }
                if (stop_insert_and_burst) {
                    hnode->burst_by_elements(curKV, cur_trie_node, hm,
                                             prefix + it->first);
                    delete hnode;
                }
            }
            return;
        }

        /*----------------searching in hash_node----------------*/

        std::pair<bool, T> find_kv_in_pages(slot* s, const CharT* key,
                                            size_t keysize) const {
            if (myTrie::hashRelative::keyEqual(get_tail_pointer(s), s->length,
                                               key, keysize)) {
                return std::pair<bool, T>(
                    true, *((T*)(get_tail_pointer(s) + s->length)));
            }
            return std::pair<bool, T>(false, T());
        }

        // return <found?, iterator>
        // iterator:    if slotid==-1, bucket is full
        //              if slotid!=-1, slotid is the insert position
        std::pair<bool, iterator> find_in_bucket(size_t bucketid,
                                                 const CharT* key,
                                                 size_t keysize) {
            // find the hitted slot in hashnode
            for (int i = 0; i != Associativity; i++) {
                slot* target_slot = get_slot(bucketid, i);
                if (target_slot->isEmpty()) {
                    return std::pair<bool, iterator>(
                        false, iterator(false, T(), this, bucketid, i));
                }

                std::pair<bool, T> res =
                    find_kv_in_pages(target_slot, key, keysize);
                if (res.first) {
                    return std::pair<bool, iterator>(
                        true, iterator(true, res.second, this, bucketid, i));
                }
            }
            return std::pair<bool, iterator>(
                false, iterator(false, T(), this, bucketid, -1));
        }

        iterator search_kv_in_hashnode(const CharT* key, size_t keysize) {
            // if found the existed target in bucket1 or bucket2, just
            // return the iterator for being modified or read
            size_t bucketId1 =
                myTrie::hashRelative::hash(key, keysize, 1) % Bucket_num;
            std::pair<bool, iterator> res1 =
                find_in_bucket(bucketId1, key, keysize);

            if (res1.first) {
                return res1.second;
            }

            size_t bucketId2 =
                myTrie::hashRelative::hash(key, keysize, 2) % Bucket_num;
            std::pair<bool, iterator> res2 =
                find_in_bucket(bucketId2, key, keysize);

            if (res2.first) {
                return res2.second;
            }

            // if the code reach here it means the target doesn't exist
            // we return the iterator with empty slot
            if (res1.second.slotid != -1) {
                return res1.second;
            } else if (res2.second.slotid != -1) {
                return res2.second;
            }
            // if two bucket are both full, we return the res1's iterator
            // with slotid == -1, and let the
            // 1. findMode: found==false
            // 2. !findMode:found==false, slotid==-1, need to kick some slot
            return res1.second;
        }

        /*----------------inserting in hash_node----------------*/

        void append_impl(const CharT* key, size_t keysize,
                         CharT* buffer_append_pos, T& value) {
            // append the string
            std::memcpy(buffer_append_pos, key, keysize * sizeof(CharT));
            buffer_append_pos += keysize;

            // append the value
            std::memcpy(buffer_append_pos, &value, sizeof(T));
        }

        // return: page_id, pos
        std::pair<size_t, size_t> alloc_insert_space(size_t keysize) {
            std::pair<char*, size_t> res = pages[cur_page_id];

            size_t need_size = keysize * sizeof(CharT) + sizeof(T);
            size_t offset = res.second;

            // if the cur_page is full, malloc a new page
            if (offset + need_size > Max_bytes_per_kv) {
                if (need_size <= Max_bytes_per_kv) {
                    char* page = (char*)malloc(Max_bytes_per_kv);
                    // set up the page information
                    pages.push_back(std::pair<char*, size_t>(page, need_size));
                    cur_page_id++;
                    return std::pair<size_t, size_t>(cur_page_id, 0);
                }

                // the need_size is surpass the max_byte_per_kv
                char* realloc_page = (char*)realloc(pages[cur_page_id].first,
                                                    offset + need_size);
                pages[cur_page_id].first = realloc_page;
                pages[cur_page_id].second = offset + need_size;

                size_t ret_page_id = cur_page_id;

                return pair<size_t, size_t>(cur_page_id, offset);
            }
            // update the page information
            pages[cur_page_id].second += need_size;
            return std::pair<size_t, size_t>(cur_page_id, offset);
        }

        std::pair<bool, T> insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, size_t bucketid,
                                                 int slotid,
                                                 std::string prefix) {
            // if slotid==-1, it denotes that the bucket(bucketid) is full ,
            // so we rehash the key_metas
            if (slotid == -1) {
                if ((slotid = cuckoo_hash(bucketid, hm)) == -1) {
                    return std::pair<bool, T>(false, T());
                }
            }

            // now the slotid cannot be -1 and slotid is lower than
            // Associativity
            assert(slotid != -1 && slotid >= 0 && slotid < Associativity);

            slot* target_slot = get_slot(bucketid, slotid);

            // allocate new page or alloc more space in old page
            std::pair<size_t, size_t> res = alloc_insert_space(keysize);

            target_slot->set_slot(keysize, res.second, res.first);
            append_impl(key, keysize, get_tail_pointer(target_slot), v);

            // set v2k
            hm->set_v2k(v, this, get_index(target_slot));
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
        anode* node;
        int index;

        SearchPoint() : node(nullptr), index(-1) {}
        SearchPoint(anode* n, int i) : node(n), index(i) {}

        void set_index(int i) { index = i; }

        std::string get_string() {
            if (node == nullptr) return string();
            // if the node is trie_node, just return the prefix on node
            if (node->is_trie_node()) {
                return ((trie_node*)node)->get_prefix();
            }

            // get the parent char chain
            trie_node* cur_node = ((hash_node*)node)->anode::parent;
            char* buf = (char*)malloc(longest_string_size);

            size_t len = cur_node->get_prefix(buf);

            // get tail
            if (index != -1) {
                class hash_node::slot* sl = ((hash_node*)node)->get_slot(index);
                memcpy(buf + len, ((hash_node*)node)->get_tail_pointer(sl),
                       sl->length);
                len += sl->length;
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
              size_t customized_bucket_count = DEFAULT_Bucket_num,
              size_t customized_byte_per_kv = DEFAULT_Max_bytes_per_kv,
              double customized_burst_ratio = DEFAULT_Burst_ratio)
        : t_root(nullptr) {
        std::cout << "SET UP GROWING-CUCKOOHASH-TRIE MAP\n";
        cout << "NO GROW\n";

        Associativity = customized_associativity;
        Bucket_num = customized_bucket_count;
        Max_bytes_per_kv = customized_byte_per_kv;
        Burst_ratio = customized_burst_ratio;

        Max_slot_num = Associativity * Bucket_num;
        Max_loop = Max_slot_num * 0.5;

        t_root = new hash_node(nullptr, string());
    }

    void set_searchPoint_index(T v, int index) { v2k[v].set_index(index); }

    void set_v2k(T v, anode* node, int index) {
        v2k[v] = SearchPoint(node, index);
    }

    void setRoot(anode* node) { t_root = node; }

    class iterator {
       public:
        bool found;
        const T v;
        anode* target_node;
        size_t bucketid;
        int slotid;

        iterator(bool f, T vv, anode* hnode, size_t bid, int sid)
            : found(f), v(vv), target_node(hnode), bucketid(bid), slotid(sid) {}

        std::pair<bool, T> insert_hashnode(const CharT* key, size_t key_size,
                                           htrie_map<CharT, T>* hm, T v) {
            return ((hash_node*)target_node)
                ->insert_kv_in_hashnode(
                    key, key_size, hm, v, bucketid, slotid,
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