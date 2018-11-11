#include <stddef.h>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include <iostream>
#include <sstream>
#include <stack>
#include "hashFunc.hpp"

#include <fstream>
#include <ios>
#include <iostream>
#include <string>

#define BURST_POINT 800
#define BUCKET_INIT_COUNT 5;

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
        static const size_t Bucket_num = 1;
        static const size_t Max_slot_num = Associativity * Bucket_num;
        static const size_t Max_bytes_per_kv = 1000;
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
            uint8_t page_id;

            bool isEmpty() { return length == 0; }
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

        iterator search_kv_in_hashnode(const CharT* key, size_t keysize) {
            if (keysize == 0) {
                if (haveValue) {
                    return iterator(true, onlyValue, this, 0, 0);
                } else {
                    return iterator(false, T(), this, 0, 0);
                }
            }
            size_t bucketId =
                myTrie::hashRelative::hash(key, keysize) % Bucket_num;

            slot* bucket_addr = (key_metas + bucketId * Associativity);
            // find the hitted slot in hashnode
            size_t i = 0;
            for (; i != Associativity; i++) {
                if (bucket_addr[i].isEmpty()) break;

                std::pair<bool, T> res =
                    find_kv_in_pages(bucket_addr[i], key, keysize);
                if (res.first) {
                    return iterator(true, res.second, this, bucketId, i);
                }
            }
            return iterator(false, T(), this, bucketId, i);
        }

        inline slot* get_slot_addr(size_t bucketid, size_t slotid) {
            return &key_metas[bucketid * Associativity + slotid];
        }

        bool need_burst() const { return elem_num >= Max_slot_num; }

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
                    target_it.insertKV(temp.data(), temp.size(), hm,
                                       itt->second);
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

        std::pair<bool, T> insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, size_t bucketid,
                                                 size_t slotid) {
            if (keysize == 0) {
                haveValue = true;
                onlyValue = v;
                hm->set_v2k(v, this, nullptr);
                return std::pair<bool, T>(true, v);
            }
            slot* target_slot = key_metas + Associativity * bucketid + slotid;

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
            if (need_burst() || slotid == Associativity - 1) {
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
        size_t slotid;

        iterator(bool f, T vv, hash_node* hnode, size_t bid, size_t sid)
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
                    return it.insertKV(key + pos, key_size - pos, this, v);
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