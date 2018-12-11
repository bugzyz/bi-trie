#include <iostream>
#include <sstream>
#include <stack>
#include "../util/hashFunc.hpp"

#include <fstream>
#include <ios>
#include <iostream>
#include <string>

// #define BUCKET_INIT_COUNT 32
// #define BURST_POINT 16384

#define BUCKET_INIT_COUNT 5
#define BURST_POINT 80

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
    class array_bucket;

    class anode {
       public:
        node_type _node_type;
        trie_node* parent;
        CharT myChar;

        bool isHashNode() { return _node_type == node_type::HASH_NODE; }
        bool isTrieNode() { return _node_type == node_type::TRIE_NODE; }
        void setParent(trie_node* p) { parent = p; }
        void setChar(CharT c) { myChar = c; }

        void deleteMe() {
            if (this->isTrieNode()) {
                std::map<CharT, trie_node*> cs = ((trie_node*)this)->childs;
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

        trie_node(CharT c, trie_node* p) {
            anode::_node_type = node_type::TRIE_NODE;
            anode::myChar = c;
            anode::parent = p;
        }

        ~trie_node() {
            std::map<CharT, trie_node*> empty;
            childs.swap(empty);
            free(onlyHashNode);
        }

        trie_node* create_trienode_With(CharT key,
                                        size_t customized_burst_threshold,
                                        size_t customized_bucket_count) {
            trie_node* new_trie_node = new trie_node(key, this);
            new_trie_node->onlyHashNode = new hash_node(
                customized_burst_threshold, customized_bucket_count);
            new_trie_node->onlyHashNode->anode::parent = new_trie_node;
            this->addChildTrieNode(new_trie_node);
            return new_trie_node;
        }

        anode* findChildNode(CharT c) {
            auto found = childs.find(c);
            if (found != childs.end()) {
                return found->second;
            } else {
                if (childs.size() == 0) {
                    return onlyHashNode;
                }
                return nullptr;
            }
        }

        void setOnlyHashNode(hash_node* node) { onlyHashNode = node; }

        hash_node* getOnlyHashNode() { return onlyHashNode; }

        void addChildTrieNode(trie_node* node) {
            childs[node->anode::myChar] = node;
            // clear the onlyHashNode because a trie_node will only have childs
            // or have the onlyHashNode
            onlyHashNode = nullptr;
        }
    };

    class hash_node : public anode {
       public:
        std::vector<array_bucket> kvs;
        T onlyValue;
        bool haveValue;
        size_t hn_burst_threshold;
        size_t hn_bucket_num;
        uint32_t element_count;

        hash_node(size_t customized_burst_threshold,
                  size_t customized_bucket_count) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = nullptr;

            hn_burst_threshold = customized_burst_threshold;
            hn_bucket_num = customized_bucket_count;

            kvs = std::vector<array_bucket>(hn_bucket_num);

            onlyValue = T();
            haveValue = false;

            element_count = 0;
        }

        ~hash_node() {
            std::vector<array_bucket> empty;
            for (auto it = kvs.begin(); it != kvs.end(); it++) {
                free(it->arr_buffer);
            }
            kvs.swap(empty);
        }

        bool need_burst() { return element_count > hn_burst_threshold; }

        T& access_onlyValue_in_hashnode(
            htrie_map* hm,
            typename myTrie::htrie_map<CharT, T>::SearchPoint* sp) {
            haveValue = true;
            if (sp != nullptr) {
                sp->hnode = this;
                sp->abucket = nullptr;
                sp->pos = 0;
            }
            return onlyValue;
        }

        // find the key with keysize in its kvs
        T& access_kv_in_hashnode(
            const CharT* key, size_t keysize, htrie_map* hm,
            typename myTrie::htrie_map<CharT, T>::SearchPoint* sp) {
            hash_node* targetNode = this;
            size_t move_pos = 0;

            // burst if hash_node have too many element
            if (need_burst()) {
                std::map<std::string, T> elements;
                // add new or existed key, if existed, it will be overwritten at
                // the loop below
                elements[std::string(key, keysize)] = T{};

                // gather elements in 'this' hash_node to burst
                for (auto it = kvs.begin(); it != kvs.end(); it++) {
                    std::vector<std::pair<std::string, T>> bucket_element;
                    it->get_item_in_array_bucket(bucket_element);
                    for (auto itt = bucket_element.begin();
                         itt != bucket_element.end(); itt++) {
                        elements[itt->first] = itt->second;
                    }
                }

                // update the targetNode
                targetNode = burst(elements, key, keysize, this->anode::parent,
                                   move_pos, hm);

                // TODO: if have onlyValue, it has to be moved to the created
                // trie_node
                delete this;

                if (targetNode->haveValue != false) {
                    if (sp != nullptr) {
                        sp->hnode = targetNode;
                        sp->abucket = nullptr;
                        sp->pos = 0;
                    }
                    return targetNode->onlyValue;
                }
            }

            size_t hashval = myTrie::hashRelative::hash<CharT>(
                                 key + move_pos, keysize - move_pos) %
                             hn_bucket_num;
            if (sp != nullptr) {
                sp->hnode = targetNode;
                return targetNode->kvs[hashval].access_kv_in_bucket(
                    key + move_pos, keysize - move_pos, sp, hm,
                    &(targetNode->element_count));
            } else {
                return targetNode->kvs[hashval].access_kv_in_bucket(
                    key + move_pos, keysize - move_pos, nullptr, hm,
                    &(targetNode->element_count));
            }
        }

        std::pair<std::string, T> getSubStr(std::pair<std::string, T>& element,
                                            size_t sub_pos = 1) {
            return std::pair<std::string, T>(element.first.substr(sub_pos),
                                             element.second);
        }

        // To turn this(a hashnode) to n childs of trie_node linking their
        // child(hashnode or trienode)
        hash_node* burst(std::map<std::string, T>& elements, const CharT* key,
                         size_t keysize, trie_node* p, size_t& move_pos,
                         htrie_map* hm) {
            hash_node* target_hashNode = nullptr;

            std::map<CharT, std::map<std::string, T>> splitElements;
            for (auto it = elements.begin(); it != elements.end(); it++) {
                splitElements[(it->first)[0]][it->first.substr(1)] = it->second;
            }

            for (auto it = splitElements.begin(); it != splitElements.end();
                 it++) {
                trie_node* cur_trie_node = new trie_node(it->first, nullptr);

                if (p == nullptr) {
                    // bursting in a root hashnode
                    // the t_root is update to a empty trie_node
                    hm->t_root = new trie_node('\0', nullptr);
                    cur_trie_node->anode::setParent((trie_node*)hm->t_root);
                    ((trie_node*)hm->t_root)->addChildTrieNode(cur_trie_node);
                    p = (trie_node*)hm->t_root;
                } else {
                    // bursting in a normal hashnode
                    cur_trie_node->anode::setParent(p);
                    p->addChildTrieNode(cur_trie_node);
                }

                if (it->first == *(key)) {
                    move_pos++;
                }

                std::map<std::string, T>& curKV = it->second;

                if (splitElements.size() == 1) {
                    return burst(curKV, key + 1, keysize - 1, cur_trie_node,
                                 move_pos, hm);
                }

                hash_node* hnode =
                    new hash_node(hn_burst_threshold, hn_bucket_num);
                hnode->anode::setParent(cur_trie_node);
                cur_trie_node->setOnlyHashNode(hnode);

                for (auto itt = curKV.begin(); itt != curKV.end(); itt++) {
                    std::string temp = itt->first;

                    if (temp.size() == 0) {
                        hnode->haveValue = true;
                        hnode->onlyValue = itt->second;
                        if (*(CharT*)key == cur_trie_node->anode::myChar &&
                            myTrie::hashRelative::keyEqual(temp.data(),
                                                           temp.size(), key + 1,
                                                           keysize - 1)) {
                            target_hashNode = hnode;
                        }
                        SearchPoint new_sp(hnode, nullptr, 0);
                        hm->v2k[itt->second] = new_sp;
                        continue;
                    }

                    size_t hashval = myTrie::hashRelative::hash<CharT>(
                                         temp.data(), temp.size()) %
                                     hn_bucket_num;
                    SearchPoint new_sp;
                    new_sp.hnode = hnode;
                    // write the value to the entry
                    hnode->kvs[hashval].access_kv_in_bucket(
                        temp.data(), temp.size(), &new_sp, hm,
                        &(hnode->element_count)) = itt->second;
                    hm->v2k[itt->second] = new_sp;

                    if (*(CharT*)key == cur_trie_node->anode::myChar &&
                        myTrie::hashRelative::keyEqual(temp.data(), temp.size(),
                                                       key + 1, keysize - 1)) {
                        target_hashNode = hnode;
                    }
                }
            }
            return target_hashNode;
        }
    };

    class array_bucket {
       public:
        CharT* arr_buffer;
        size_t buffer_size;
        static const KeySizeT END_OF_BUCKET =
            std::numeric_limits<KeySizeT>::max();

        array_bucket() {
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
                if (myTrie::hashRelative::keyEqual(cmp_buffer_ptr, length, key,
                                                   keysize)) {
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

        void resetBucketElement_v2k(CharT* arr_buffer,
                                    htrie_map<CharT, T>* hm) {
            CharT* buffer_ptr = arr_buffer;
            size_t pos = 0;
            while (!is_end_of_bucket(buffer_ptr)) {
                size_t length = read_key_size(buffer_ptr);
                T* value_pos = (T*)(buffer_ptr + sizeof(KeySizeT) + length);

                (hm->v2k[*value_pos]).pos = buffer_ptr;

                // move ptr to next header, skip keysize, string, value
                buffer_ptr = buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
            }
        }

        T& access_kv_in_bucket(const CharT* target, size_t keysize,
                               typename htrie_map<CharT, T>::SearchPoint* sp,
                               htrie_map<CharT, T>* hm, uint32_t* counter) {
            std::pair<size_t, bool> found = find_in_bucket(target, keysize);

            if (found.second) {
                // found the target, return the value reference
                if (sp != nullptr) {
                    sp->abucket = this;
                    sp->pos = found.first;
                }
                return get_val_ref(found.first);

            } else {
                // not found, the buffer is full, need realloc
                size_t new_size = cal_arrbuffer_newsize(keysize);

                CharT* new_buffer =
                    (CharT*)(std::realloc(arr_buffer, new_size));
                if (new_buffer == nullptr) {
                    std::cerr << "realloc failed!!!!!!!!!1" << std::endl;
                    exit(-1);
                }
                arr_buffer = new_buffer;
                buffer_size = new_size;

                T v = T{};
                append_impl(target, keysize, arr_buffer + found.first, v);
                if (sp != nullptr) {
                    sp->abucket = this;
                    sp->pos = found.first;
                }
                (*counter)++;
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

        void get_item_in_array_bucket(
            std::vector<std::pair<std::string, T>>& res) {
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
            }
        }
    };

    class SearchPoint {
       public:
        hash_node* hnode;
        array_bucket* abucket;
        uint32_t pos;

        SearchPoint() : hnode(nullptr), pos(0) {}
        SearchPoint(hash_node* h, array_bucket* ab, uint32_t p)
            : hnode(h), abucket(ab), pos(p) {}

        std::string getString() {
            trie_node* cur_node = hnode->anode::parent;
            std::stack<char> st;
            while (cur_node != nullptr && cur_node->anode::myChar != '\0') {
                st.push((char)cur_node->anode::myChar);
                cur_node = cur_node->anode::parent;
            }
            std::stringstream ss;
            while (!st.empty()) {
                ss << st.top();
                st.pop();
            }
            if (abucket != nullptr) {
                std::string res;
                KeySizeT key_size;
                CharT* buf = abucket->arr_buffer;
                std::memcpy(&key_size, buf + pos, sizeof(KeySizeT));

                size_t length = (size_t)key_size;

                char* temp = (char*)malloc(length + 1);
                std::memcpy(temp, buf + pos + sizeof(KeySizeT), length);
                temp[length] = '\0';
                res = std::string(temp);
                free(temp);
                ss << res;
            }
            return ss.str();
        }
    };

   public:
    anode* t_root;
    std::map<T, SearchPoint> v2k;
    htrie_map(size_t customized_burst_threshold = DEFAULT_BURST_THRESHOLD,
              size_t customized_bucket_count = DEFAULT_BURST_THRESHOLD) {
        std::cout << "SET UP HAT-TRIE MAP\n";
        t_root =
            new hash_node(customized_burst_threshold, customized_bucket_count);
        burst_threshold = customized_burst_threshold;
        bucket_num = customized_bucket_count;
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
            return std::pair<bool, T>(false, std::string());
        } else {
            return std::pair<bool, T>(true, v2k[v].get_string());
        }
    }

    /*----------------------------------*/

    bool insertKV(std::string key, T v) {
        SearchPoint sp;
        find(key.data(), key.size(), &sp) = v;
        v2k[v] = sp;
        return true;
    }

    T& find(const CharT* key, size_t key_size, SearchPoint* sp) {
        anode* current_node = t_root;
        if (t_root->isHashNode()) {
            // find a key in hash_node
            // if find: return the v reference
            // if notfind: write the key and return the v reference
            return ((hash_node*)current_node)
                ->access_kv_in_hashnode(key, key_size, this, sp);
        }

        for (size_t pos = 0; pos < key_size; pos++) {
            if (current_node->isTrieNode()) {
                anode* parent = current_node;
                current_node =
                    ((trie_node*)current_node)->findChildNode(key[pos]);

                if (current_node == nullptr) {
                    // can't find, create a relative trie_node and a hashnode
                    // child, set current_node as the trie_node
                    current_node = ((trie_node*)parent)
                                       ->create_trienode_With(
                                           key[pos], this->burst_threshold,
                                           this->bucket_num);
                } else if (current_node->isHashNode()) {
                    // if find the target hashnode instead of moving the pos,
                    // pos should be recover
                    pos--;
                }
            } else {
                return ((hash_node*)current_node)
                    ->access_kv_in_hashnode(key + pos, key_size - pos, this,
                                            sp);
            }
        }
        return ((trie_node*)current_node)
            ->getOnlyHashNode()
            ->access_onlyValue_in_hashnode(this, sp);
    }

    void deleteMyself() { t_root->deleteMe(); }
};  // namespace myTrie

}  // namespace myTrie