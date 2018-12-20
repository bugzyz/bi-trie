#include <stddef.h>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include <iostream>
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

// todo: wait to be deleted, just for recording the time that expand() cost
uint64_t rehash_cost_time = 0;
uint64_t rehash_total_num = 0;

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
    enum class node_type : unsigned char { HASH_NODE, TRIE_NODE };
    class trie_node;
    class hash_node;
    class iterator;

    class anode {
       public:
        node_type _node_type;
        trie_node* parent;

        bool is_hash_node() { return _node_type == node_type::HASH_NODE; }
        bool is_trie_node() { return _node_type == node_type::TRIE_NODE; }
        void setParent(trie_node* p) { parent = p; }

        void deleteMe() {
            if (this->is_trie_node()) {
                std::map<CharT, trie_node*> cs = ((trie_node*)this)->childs;
                for (auto it = cs.begin(); it != cs.end(); it++) {
                    it->second->deleteMe();
                    delete it->second;
                }
            } else {
                delete (hash_node*)this;
            }
        }
    };

    class trie_node : public anode {
       public:
        // store the suffix of hash_node or trie_node
        std::map<CharT, trie_node*> childs;
        hash_node* onlyHashNode;
        CharT myChar;

        trie_node(CharT c, trie_node* p) {
            anode::_node_type = node_type::TRIE_NODE;
            anode::parent = p;

            myChar = c;
            onlyHashNode = nullptr;
        }

        ~trie_node() {
            std::map<CharT, trie_node*> empty;
            childs.swap(empty);
        }

        anode* findChildNode(CharT c, bool findMode) {
            auto found = childs.find(c);
            if (found != childs.end()) {
                trie_node* target = found->second;
                if (target->onlyHashNode != nullptr)
                    return target->onlyHashNode;
                return found->second;
            } else {
                if (findMode) {
                    return nullptr;
                } else {
                    trie_node* son_trie_node = new trie_node(c, this);
                    this->addChildTrieNode(son_trie_node);
                    son_trie_node->onlyHashNode = new hash_node(son_trie_node);
                    return son_trie_node->onlyHashNode;
                }
            }
        }

        void setOnlyHashNode(hash_node* node) { onlyHashNode = node; }

        hash_node* getOnlyHashNode() { return onlyHashNode; }

        void addChildTrieNode(trie_node* node) {
            childs[node->myChar] = node;
            // clear the onlyHashNode because a trie_node will only have childs
            // or have the onlyHashNode
            onlyHashNode = nullptr;
        }
    };

    class hash_node : public anode {
       public:
        class slot;

        slot* key_metas;
        vector<std::pair<char*, size_t>> pages;
        size_t elem_num;
        size_t cur_page_id;

        T onlyValue;
        bool haveValue;

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
        };

       public:
        explicit hash_node(trie_node* p)
            : elem_num(0), cur_page_id(0), onlyValue(T()), haveValue(false) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = p;
            key_metas = (slot*)malloc(Max_slot_num * sizeof(slot));

            // init key space
            for (int i = 0; i != Max_slot_num; i++) {
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

        inline char* get_tail_pointer(slot* s) {
            return pages[s->page_id].first + s->pos;
        }

        std::pair<bool, T> find_kv_in_pages(slot& s, const CharT* key,
                                            size_t keysize) const {
            if (myTrie::hashRelative::keyEqual(pages[s.page_id].first + s.pos,
                                               s.length, key, keysize)) {
                return std::pair<bool, T>(
                    true, *((T*)(pages[s.page_id].first + s.pos + s.length)));
            }
            return std::pair<bool, T>(false, T());
        }

        // return <found?, iterator>
        // iterator:    if slotid==-1, bucket is full
        //              if slotid!=-1, slotid is the insert position
        std::pair<bool, iterator> find_in_bucket(size_t bucketid,
                                                 const CharT* key,
                                                 size_t keysize) {
            slot* bucket_addr = (key_metas + bucketid * Associativity);
            // find the hitted slot in hashnode
            for (int i = 0; i != Associativity; i++) {
                if (bucket_addr[i].isEmpty()) {
                    return std::pair<bool, iterator>(
                        false, iterator(false, T(), this, bucketid, i));
                }

                std::pair<bool, T> res =
                    find_kv_in_pages(bucket_addr[i], key, keysize);
                if (res.first) {
                    return std::pair<bool, iterator>(
                        true, iterator(true, res.second, this, bucketid, i));
                }
            }
            return std::pair<bool, iterator>(
                false, iterator(false, T(), this, bucketid, -1));
        }

        iterator search_kv_in_hashnode(const CharT* key, size_t keysize) {
            if (keysize == 0) {
                return iterator(haveValue, onlyValue, this, 0, 0);
            }

            size_t bucketId1 =
                myTrie::hashRelative::hash(key, keysize, 1) % Bucket_num;
            std::pair<bool, iterator> res1 =
                find_in_bucket(bucketId1, key, keysize);

            size_t bucketId2 =
                myTrie::hashRelative::hash(key, keysize, 2) % Bucket_num;
            std::pair<bool, iterator> res2 =
                find_in_bucket(bucketId2, key, keysize);

            // if found the existed target in bucket1 or bucket2, just return
            // the iterator for being modified or read
            if (res1.first) {
                return res1.second;
            } else if (res2.first) {
                return res2.second;
            }

            // if the code reach here it means the target doesn't exist
            // we return the iterator with empty slot
            if (res1.second.slotid != -1) {
                return res1.second;
            } else if (res2.second.slotid != -1) {
                return res2.second;
            }
            // if two bucket are both full, we return the res1's iterator with
            // slotid == -1, and let the
            // 1. findMode: found==false
            // 2. !findMode:found==false, slotid==-1, need to kick some slot
            return res1.second;
        }

        inline slot* get_slot_addr(size_t bucketid, size_t slotid) {
            return &key_metas[bucketid * Associativity + slotid];
        }

        bool need_burst() const {
            return elem_num >= Max_slot_num * Burst_ratio;
        }

        // To turn this(a hashnode) to n childs of trie_node linking their
        // hashnode
        void burst(std::map<std::string, T>& elements, trie_node* p,
                   htrie_map* hm) {
            std::map<CharT, std::map<std::string, T>> preprocElements;
            for (auto it = elements.begin(); it != elements.end(); it++) {
                preprocElements[(it->first)[0]][it->first.substr(1)] =
                    it->second;
            }

            for (auto it = preprocElements.begin(); it != preprocElements.end();
                 it++) {
                if (p == nullptr) {
                    // bursting in a root hashnode
                    // the t_root is update to a empty trie_node
                    p = new trie_node('\0', nullptr);
                    hm->setRoot(p);
                }

                trie_node* cur_trie_node = new trie_node(it->first, p);
                p->addChildTrieNode(cur_trie_node);

                std::map<std::string, T>& curKV = it->second;
                if (preprocElements.size() == 1) {
                    return burst(curKV, cur_trie_node, hm);
                }

                hash_node* hnode = new hash_node(cur_trie_node);
                cur_trie_node->setOnlyHashNode(hnode);

                bool stop_insert_and_burst = false;
                for (auto itt = curKV.begin(); itt != curKV.end(); itt++) {
                    std::string temp = itt->first;

                    if (temp.size() == 0) {
                        hnode->haveValue = true;
                        hnode->onlyValue = itt->second;

                        hm->set_v2k(itt->second, hnode, nullptr);
                        continue;
                    }

                    iterator target_it =
                        hnode->search_kv_in_hashnode(temp.data(), temp.size());
                    std::pair<bool, T> res = target_it.insertKV(
                        temp.data(), temp.size(), hm, itt->second);
                    // if insert failed, it need burst
                    if (res.first == false) {
                        stop_insert_and_burst = true;
                        break;
                    }
                }
                if (stop_insert_and_burst) {
                    burst(curKV, cur_trie_node, hm);
                    delete hnode;
                }
            }
            return;
        }

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

            // if the cur_page is full, malloc a new page
            if (res.second + keysize * sizeof(CharT) + sizeof(T) >
                Max_bytes_per_kv) {
                char* page = (char*)malloc(Max_bytes_per_kv);
                // set up the page information
                pages.push_back(std::pair<char*, size_t>(
                    page, keysize * sizeof(CharT) + sizeof(T)));
                cur_page_id++;
                return std::pair<size_t, size_t>(cur_page_id, 0);
            }
            size_t offset = res.second;
            // update the page information
            pages[cur_page_id].second += keysize * sizeof(CharT) + sizeof(T);
            return std::pair<size_t, size_t>(cur_page_id, offset);
        }

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

        inline slot* get_bucket(size_t bucketid) {
            return key_metas + bucketid * Associativity;
        }

        inline slot* get_slot(size_t bucketid, size_t slotid) {
            return key_metas + bucketid * Associativity + slotid;
        }

        inline slot* previous_dst_slot_in_same_bucket(slot* s) {
            size_t slotid = (s - key_metas) % Associativity;
            if (slotid == 0) {
                return nullptr;
            } else
                return s - 1;
        }

        void apply_the_changed_searchPoint(map<T, slot*>& searchPoints,
                                           htrie_map<CharT, T>* hm) {
            for (auto it = searchPoints.begin(); it != searchPoints.end(); it++)
                hm->set_searchPoint_slot(it->first, it->second);
        }

        int rehash(size_t bucketid, htrie_map<CharT, T>* hm) {
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

            map<T, slot*> searchPoint_wait_2_be_update;
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
                        dst_slot;

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
                        dst_slot;
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

        std::pair<bool, T> insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, size_t bucketid,
                                                 int slotid) {
            if (keysize == 0) {
                haveValue = true;
                onlyValue = v;
                hm->set_v2k(v, this, nullptr);
                return std::pair<bool, T>(true, v);
            }

            // if slotid==-1, it denotes that the bucket(bucketid) is full , so
            // we rehash the key_metas
            if (slotid == -1) {
                if ((slotid = rehash(bucketid, hm)) == -1) {
                    return std::pair<bool, T>(false, T());
                }
            }

            // now the slotid cannot be -1 and slotid is lower than
            // Associativity
            assert(slotid != -1 && slotid >= 0 && slotid < Associativity);

            slot* target_slot =
                key_metas + Associativity * bucketid + (size_t)slotid;

            // allocate new page or alloc more space in old page
            std::pair<size_t, size_t> res = alloc_insert_space(keysize);

            target_slot->length = keysize;
            target_slot->page_id = res.first;
            target_slot->pos = res.second;
            append_impl(key, keysize, get_tail_pointer(target_slot), v);

            // set v2k
            hm->set_v2k(v, this, target_slot);
            elem_num++;

            // todo: need to burst elegantly
            if (need_burst()) {
                std::map<std::string, T> elements;
                get_all_elements(elements);

                burst(elements, this->anode::parent, hm);

                delete this;
            }

            return std::pair<bool, T>(true, v);
        }
    };

    class SearchPoint {
       public:
        hash_node* hnode;
        typename hash_node::slot* sl;

        SearchPoint() : hnode(nullptr), sl(nullptr) {}
        SearchPoint(hash_node* h, typename hash_node::slot* s)
            : hnode(h), sl(s) {}

        void set_slot(typename hash_node::slot* s) { sl = s; }

        std::string get_string() {
            // get the parent char chain
            trie_node* cur_node = hnode->anode::parent;
            std::string res;
            while (cur_node != nullptr && cur_node->myChar != '\0') {
                res = (char)cur_node->myChar + res;
                cur_node = cur_node->anode::parent;
            }
            // get tail
            if (sl != nullptr) {
                res = res + std::string(hnode->get_tail_pointer(sl),
                                        (size_t)sl->length);
            }
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
        std::cout << "SET UP CUCKOOHASH-TRIE MAP\n";

        Associativity = customized_associativity;
        Bucket_num = customized_bucket_count;
        Max_bytes_per_kv = customized_byte_per_kv;
        Burst_ratio = customized_burst_ratio;

        Max_slot_num = Associativity * Bucket_num;
        Max_loop = Max_slot_num * 0.5;

        t_root = new hash_node(nullptr);
    }

    void set_searchPoint_slot(T v, typename hash_node::slot* s) {
        v2k[v].set_slot(s);
    }

    void set_v2k(T v, hash_node* hnode, typename hash_node::slot* s) {
        v2k[v] = SearchPoint(hnode, s);
    }

    void setRoot(anode* node) { t_root = node; }

    /*----------------------------------*/

    // access element
    T searchByKey(std::string key) {
        return access_kv_in_htrie_map(key.data(), key.size(), T(), true).second;
    }

    std::string searchByValue(T v) { return v2k[v].get_string(); }

    // find operation
    std::pair<bool, std::string> findByKey(std::string key) {
        return access_kv_in_htrie_map(key.data(), key.size(), T(), true);
    }

    std::pair<bool, T> findByValue(T v) {
        if (v2k.find(v) == v2k.end()) {
            return std::pair<bool, T>(false, string());
        } else {
            return std::pair<bool, T>(true, v2k[v].get_string());
        }
    }

    /*----------------------------------*/

    std::pair<bool, T> insertKV(std::string key, T v) {
        return access_kv_in_htrie_map(key.data(), key.size(), v, false);
    }

    class iterator {
       public:
        bool found;
        const T v;
        hash_node* target_node;
        size_t bucketid;
        int slotid;

        iterator(bool f, T vv, hash_node* hnode, size_t bid, int sid)
            : found(f), v(vv), target_node(hnode), bucketid(bid), slotid(sid) {}

        std::pair<bool, T> insertKV(const CharT* key, size_t key_size,
                                    htrie_map<CharT, T>* hm, T v) {
            return target_node->insert_kv_in_hashnode(key, key_size, hm, v,
                                                      bucketid, slotid);
        }
    };

    std::pair<bool, T> access_kv_in_htrie_map(const CharT* key, size_t key_size,
                                              T v, bool findMode) {
        anode* current_node = t_root;

        for (size_t pos = 0; pos < key_size; pos++) {
            if (current_node->is_trie_node()) {
                trie_node* parent;
                parent = (trie_node*)current_node;

                // only return the hitted trie_node* or nullptr if not found
                current_node = parent->findChildNode(key[pos], findMode);

            } else {
                iterator it =
                    ((hash_node*)current_node)
                        ->search_kv_in_hashnode(key + pos, key_size - pos);
                if (findMode) {
                    return std::pair<bool, T>(it.found, it.v);
                } else {
                    pair<bool, T> res =
                        it.insertKV(key + pos, key_size - pos, this, v);
                    if (res.first == false) {
                        // if the insert failed, we burst the target_hashnode
                        // and retry insertion
                        hash_node* hnode_burst_needed =
                            (hash_node*)current_node;
                        map<string, T> hnode_elems;
                        hnode_burst_needed->get_all_elements(hnode_elems);
                        hnode_burst_needed->burst(
                            hnode_elems, hnode_burst_needed->anode::parent,
                            this);
                        delete hnode_burst_needed;
                        return access_kv_in_htrie_map(key, key_size, v, false);
                    }
                    return res;
                }
            }
            // only in the findMode==true can cause the current_node to be
            // nullptr
            if (current_node == nullptr) {
                return std::pair<bool, T>(false, T());
            }
        }

        // find a key in hash_node's only value
        iterator it = ((hash_node*)current_node)->search_kv_in_hashnode(key, 0);

        if (findMode) {
            return std::pair<bool, T>(it.found, it.v);
        } else {
            return it.insertKV(key, 0, this, v);
        }
    }

    void deleteMyself() {
        map<T, SearchPoint> empty;
        v2k.swap(empty);
        t_root->deleteMe();
    }
};  // namespace myTrie

}  // namespace myTrie