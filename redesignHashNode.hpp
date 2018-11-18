#include <stddef.h>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include <iostream>
#include <sstream>
#include <stack>
#include "_hashFunc.hpp"

#include <fstream>
#include <ios>
#include <iostream>
#include <set>
#include <string>

#include <time.h>

#define BURST_POINT 800
#define BUCKET_INIT_COUNT 5;

// debug
int rehashFailed = 0;
int failed_at_the_first_place = 0;

using namespace std;

namespace myTrie {
// charT = char, T = value type, keysizeT = the type describe keysize
template <class CharT, class T, class KeySizeT = std::uint16_t>
class htrie_map {
   public:
    // burst_threshold should be set to be greater than 26 for 26 alaphbet
    // and ", @ .etc.(the test in lubm40 shows that it has 50 char species)
    static const size_t DEFAULT_BURST_THRESHOLD = BURST_POINT;
    static const size_t DEFAULT_BUCKET_INIT_COUNT = BUCKET_INIT_COUNT;
    size_t burst_threshold;
    size_t bucket_num;

    enum class node_type : unsigned char { HASH_NODE, TRIE_NODE };
    class trie_node;
    class hash_node;
    class iterator;

    class anode {
       public:
        node_type _node_type;
        trie_node* parent;

        bool isHashNode() { return _node_type == node_type::HASH_NODE; }
        bool isTrieNode() { return _node_type == node_type::TRIE_NODE; }
        void setParent(trie_node* p) { parent = p; }

        void deleteMe() {
            if (this->isTrieNode()) {
                std::map<CharT, anode*> cs = ((trie_node*)this)->childs;
                for (auto it = cs.begin(); it != cs.end(); it++) {
                    it->second->deleteMe();
                    delete it->second;
                }
            } else {
                delete this;
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
            std::map<CharT, anode*> empty;
            childs.swap(empty);
            free(onlyHashNode);
        }

        hash_node* create_trienode_With(CharT key,
                                        size_t customized_burst_threshold,
                                        size_t customized_bucket_count) {
            trie_node* new_trie_node = new trie_node(key, this);
            this->addChildTrieNode(new_trie_node);
            new_trie_node->onlyHashNode = new hash_node(new_trie_node);
            return new_trie_node->onlyHashNode;
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
        static const size_t Associativity = 4;
        static const size_t Bucket_num = 11;
        static const size_t Max_slot_num = Associativity * Bucket_num;
        static const size_t Max_bytes_per_kv = 1000;
        static const size_t Max_loop = Max_slot_num * 0.5;
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
            size_t pos;
            size_t page_id;

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
            size_t slot_num = Bucket_num * Associativity;
            key_metas = (slot*)malloc(slot_num * sizeof(slot));

            // init key space
            for (int i = 0; i != slot_num; i++) {
                key_metas[i].length = 0;
                key_metas[i].pos = 0;
                key_metas[i].page_id = 0;
            }

            char* page = (char*)malloc(Max_bytes_per_kv);
            pages.push_back(std::pair<char*, size_t>(page, 0));
        }

        // todo
        ~hash_node() {}

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
            int i = 0;
            for (; i != Associativity; i++) {
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
                if (haveValue) {
                    return iterator(true, onlyValue, this, 0, 0);
                } else {
                    return iterator(false, T(), this, 0, 0);
                }
            }
            size_t bucketId1 =
                myTrie::hashRelative::hash(key, keysize, 1) % Bucket_num;
            std::pair<bool, iterator> res1 =
                find_in_bucket(bucketId1, key, keysize);

            size_t bucketId2 =
                myTrie::hashRelative::hash(key, keysize, 2) % Bucket_num;
            std::pair<bool, iterator> res2 =
                find_in_bucket(bucketId1, key, keysize);

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

        bool need_burst() const { return elem_num >= Max_slot_num * 0.5; }

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
                    if (itt->second == 894180) {
                        cout << "'debug2'" << endl;
                    }
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
            // T v;
            std::memcpy(&elements[res], tail_pointer + s->length, sizeof(T));
            // elements[res] = v;
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

        int rehash(size_t bucketid) {
            cout << "========new rehash==============\n";
            print_key_metas();
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

            size_t kicked_slot_id = ret_slot_id;
            size_t current_bucket_id = bucketid;

            slot src_slot = slot(0, 0, 0);
            slot* dst_slot = get_slot(current_bucket_id, kicked_slot_id);

            size_t rehash_count = 0;
            // kicking slot
            do {
                cout << "rehash time: " << rehash_count << "\n";
                cout << "src_slot: " << src_slot.length << "/" << src_slot.pos
                     << "/" << src_slot.page_id << "\n";
                cout << "dst_num: " << dst_slot - key_metas << ": "
                     << dst_slot->length << "/" << dst_slot->pos << "/"
                     << dst_slot->page_id << endl;
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
                if (bucketid_kick_to == current_bucket_id) {
                    cout << "num." << dst_slot - key_metas
                         << " is stable in bucket: " << bucketid_kick_to
                         << endl;
                    dst_slot = dst_slot - 1;
                    ret_slot_id = ret_slot_id - 1;
                    if (ret_slot_id == -1) {
                        return -1;
                    }
                    continue;
                }

                cout << "current_bucket_id: " << current_bucket_id << "\n";
                cout << "bucketid_kick_to: " << bucketid_kick_to << "\n";

                /*
                    src(a,b,c)
                    temp: d,e,f
                    cur_bucket: |x      |x      |x      |dst(d,e,f)|
                    kk2_bucket: |x      |x      |x      |x         |
                */
                KeySizeT temp_length = dst_slot->length;
                size_t temp_pos = dst_slot->pos;
                size_t temp_page_id = dst_slot->page_id;
                cout << "temp: " << temp_length << "/" << temp_pos << "/"
                     << temp_page_id << "/"
                     << "\n";
                cout << "--------------------------------\n";

                /*
                    src(a,b,c)
                    temp: d,e,f
                    cur_bucket: |x      |x      |x      |dst(a,b,c)|
                    kk2_bucket: |x      |x      |x      |x         |
                */
                dst_slot->set_slot(src_slot.length, src_slot.pos,
                                   src_slot.page_id);

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

                    print_key_metas();
                    cout << "=======rehash finished!=========\n";

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
                     cur_bucket: |x      |x      |x      |src(j,k,l)|
                     kk2_bucket: |x      |x      |x      |x-dst     |
                */
                src_slot.set_slot(temp_length, temp_pos, temp_page_id);
                /*
                     src(d,e,f)
                     (kk2_bucket)
                     cur_bucket: |x      |x      |x      |x-dst     |
                */
                current_bucket_id = bucketid_kick_to;

                // print_key_metas();

                rehash_count++;
            } while (rehash_count != Max_loop);
            return -1;
        }

        struct SlotCmp {
            bool operator()(const slot& lhs, const slot& rhs) const {
                return lhs.length * 100 + lhs.pos * 10 + lhs.page_id <
                       rhs.length * 100 + rhs.pos * 10 + rhs.page_id;
            }
        };

        // just element
        void print_hashnode_element() {
            cout << "printing hashnode elememt\n";
            std::map<string, T> tempMap;
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != Associativity; j++) {
                    slot* sss = get_slot(i, j);
                    if (!sss->isEmpty()) {
                        get_tail_str_v(tempMap, sss);
                    }
                }
            }
            cout << "hashnode elememt_num:" << elem_num
                 << " acutally:" << tempMap.size() << "\n";
            if (elem_num != tempMap.size()) {
                assert(true);
            }
            int count = 0;
            for (auto it = tempMap.begin(); it != tempMap.end(); it++) {
                cout << count++ << ": " << it->first << " = " << it->second
                     << endl;
            }
        }

        // just key_metas layout
        void print_key_metas() {
            cout << "print keymetas layout\n";
            for (int i = 0; i != Bucket_num; i++) {
                cout << i << ":\t";
                for (int j = 0; j != Associativity; j++) {
                    slot* s = get_slot(i, j);
                    cout << i * Associativity + j << ":" << s->length << "/"
                         << s->pos << "/" << s->page_id << "\t\t";
                }
                cout << "\n";
            }
        }

        void print_slot(slot s) {
            cout << s.length << "/" << s.pos << "/" << s.page_id << endl;
        }

        void setup_before_slot_situation(set<slot, SlotCmp>& checkingset) {
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != Associativity; j++) {
                    slot* s = get_slot(i, j);
                    if (!s->isEmpty()) checkingset.insert(slot(*s));
                }
            }
        }

        void check_current_slot_situation(set<slot, SlotCmp>& checkingSet) {
            set<slot, SlotCmp> secondCheck;
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != Associativity; j++) {
                    slot* s = get_slot(i, j);
                    if (!s->isEmpty()) secondCheck.insert(slot(*s));
                }
            }
            cout << checkingSet.size() << "-" << secondCheck.size() << endl;
            if (checkingSet.size() != secondCheck.size()) {
                for (auto it = secondCheck.begin(); it != secondCheck.end();
                     it++) {
                    auto itt = checkingSet.find(*it);

                    if (itt == checkingSet.end()) {
                        cout << "cannt find: ";
                        print_slot(*it);
                    }
                }
                print_hashnode_element();
                print_key_metas();
            }
        }

        std::pair<bool, T> insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, size_t bucketid,
                                                 int slotid) {
            static hash_node* debug_hnode;
            if (v == 894180) {
                cout << "'debug'" << endl;
                debug_hnode = this;
            }

            if (keysize == 0) {
                haveValue = true;
                onlyValue = v;
                hm->set_v2k(v, this, nullptr);
                return std::pair<bool, T>(true, v);
            }

            if (this == debug_hnode) {
                cout << "check it out\n";
            }

            // if slotid==-1, it denotes that the bucket(bucketid) is full , so
            // we rehash the key_metas
            if (slotid == -1) {
                set<slot, SlotCmp> checkingSet;
                setup_before_slot_situation(checkingSet);

                if ((slotid = rehash(bucketid)) == -1) {
                    rehashFailed++;
                    check_current_slot_situation(checkingSet);

                    cout << "rehash failed!\n";
                    return std::pair<bool, T>(false, T());
                }
                cout << "Rehashing success: slotid is updated to " << slotid
                     << "\n";
                check_current_slot_situation(checkingSet);
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

            if (this == debug_hnode) {
                print_hashnode_element();
                print_key_metas();
            }

            // todo: need to burst elegantly
            if (need_burst()) {
                std::map<std::string, T> elements;
                get_all_elements(elements);

                if (this == debug_hnode) {
                    print_hashnode_element();
                    print_key_metas();
                }

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

        // todo: refine
        std::string getString() {
            // get the parent char chain
            trie_node* cur_node = hnode->anode::parent;
            std::stack<char> st;
            while (cur_node != nullptr && cur_node->myChar != '\0') {
                st.push((char)cur_node->myChar);
                cur_node = cur_node->anode::parent;
            }
            std::stringstream ss;
            while (!st.empty()) {
                ss << st.top();
                st.pop();
            }
            // get tail
            if (sl != nullptr) {
                string res = std::string(hnode->get_tail_pointer(sl),
                                         (size_t)sl->length);
                ss << res;
            }
            return ss.str();
        }
    };

   public:
    anode* t_root;
    std::map<T, SearchPoint> v2k;
    htrie_map(size_t customized_burst_threshold = DEFAULT_BURST_THRESHOLD,
              size_t customized_bucket_count = DEFAULT_BURST_THRESHOLD)
        : t_root(new hash_node(nullptr)),
          burst_threshold(customized_burst_threshold),
          bucket_num(customized_bucket_count) {}

    void set_v2k(T v, hash_node* hnode, typename hash_node::slot* s) {
        v2k[v] = SearchPoint(hnode, s);
    }

    void setRoot(anode* node) { t_root = node; }

    std::pair<bool, T> searchByKey(std::string key) {
        return access_kv_in_htrie_map(key.data(), key.size(), T(), true);
    }

    std::string searchByValue(T v) { return v2k[v].getString(); }

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
            if (current_node->isTrieNode()) {
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
                        failed_at_the_first_place++;
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

    void deleteMyself() { t_root->deleteMe(); }
};  // namespace myTrie

}  // namespace myTrie