#include <cstdint>
// for std::strlen
#include <cstring>
#include <map>
#include <vector>

#include <string>

// numeric_limit
#include <limits>

// debug
#include <fstream>
using std::fstream;
#include <iostream>

#define debugStream \
    if (true) {     \
    } else          \
        std::cout

// #define BUCKET_INIT_COUNT 32
// #define BURST_POINT 16384

#define BUCKET_INIT_COUNT 5
#define BURST_POINT 80

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
    if (key_size_lhs == 0 && key_size_rhs == 0) {
        return true;
    }
    if (key_size_lhs != key_size_rhs) {
        debugStream << "size not equal" << std::endl;
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
    // burst_threshold should be set to be greater than 26 for 26 alaphbet
    // and ", @ .etc.(the test in lubm40 shows that it has 50 char species)
    static const size_t DEFAULT_BURST_THRESHOLD = BURST_POINT;
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
        void setChar(CharT c) { myChar = c; }
    };
    class trie_node : public anode {
       public:
        // store the suffix of hash_node or trie_node
        std::map<CharT, anode*> childs;
        hash_node* onlyHashNode;

        trie_node(CharT c, trie_node* p) {
            anode::_node_type = node_type::TRIE_NODE;
            anode::myChar = c;
            anode::parent = p;
        }

        trie_node* create_trienode_With(CharT key) {
            trie_node* new_trie_node = new trie_node(key, this);
            // TODO: remove the hard code
            new_trie_node->onlyHashNode = new hash_node(BURST_POINT);
            new_trie_node->onlyHashNode->anode::parent = new_trie_node;
            this->addChildTrieNode(new_trie_node);
            debugStream << "father: " << (void*)this
                        << " create a new trie_node: " << (void*)new_trie_node
                        << " return the onlyHahsNode: "
                        << (void*)new_trie_node->onlyHashNode << std::endl;
            return new_trie_node;
        }

        std::map<CharT, anode*> getChildsMap() { return childs; }

        anode* findChildNode(CharT c) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                debugStream << (void*)this << " have child: " << it->first
                            << " at " << (void*)it->second << std::endl;
            }
            auto found = childs.find(c);
            if (found != childs.end()) {
                return found->second;
                debugStream << "return target\n";

            } else {
                if (childs.size() == 0) {
                    debugStream << "no match child return onlyHashNode: "
                                << (void*)onlyHashNode << std::endl;
                    return onlyHashNode;
                }
                debugStream << "return nullptr\n";
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

        void addChildsNode(
            std::vector<std::pair<CharT, trie_node*>> newChilds) {
            for (int i = 0; i != newChilds.size(); i++) {
                addChildsNode(newChilds[i]);
            }
        }
    };

    class hash_node : public anode {
       public:
        std::vector<array_bucket> kvs;
        size_t burst_threshold;
        T onlyValue;
        bool haveValue;

        hash_node(size_t customized_burst_threshold) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = nullptr;
            kvs = std::vector<array_bucket>(BUCKET_INIT_COUNT);
            burst_threshold = customized_burst_threshold;
            debugStream << "new hash_node is init, burst_threshold set to be "
                        << burst_threshold << std::endl;
            onlyValue = T();
            haveValue = false;
        }

        bool need_burst() {
            size_t node_element_count = 0;
            for (auto it = kvs.begin(); it != kvs.end(); it++) {
                node_element_count += it->bucket_element_count;
            }
            debugStream << "node element count: " << node_element_count
                        << std::endl;
            debugStream << "burst_threshold: " << burst_threshold << std::endl;
            return node_element_count > burst_threshold;
        }

        T& access_onlyValue_in_hashnode(htrie_map* hm) {
            haveValue = true;
            return onlyValue;
        }

        // find the key with keysize in its kvs
        T& access_kv_in_hashnode(const CharT* key, size_t keysize,
                                 htrie_map* hm) {
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

                for (auto it = elements.begin(); it != elements.end(); it++) {
                    debugStream << "elements: " << it->first
                                << " with value: " << it->second << std::endl;
                }
                debugStream << "bursting in " << (void*)this << std::endl;

                // update the targetNode
                targetNode = burst(elements, key, keysize, this->anode::parent,
                                   move_pos, hm);
                debugStream << "target_node found: " << (void*)targetNode
                            << std::endl;

                // TODO: if have onlyValue, it has to be moved to the created
                // trie_node
                delete this;

                if (targetNode->haveValue != false) {
                    return targetNode->onlyValue;
                }
            }

            size_t hashval = myTrie::hashRelative::hash<CharT>(
                key + move_pos, keysize - move_pos);
            debugStream << "bucket_id: " << hashval << std::endl;
            return targetNode->kvs[hashval].access_kv_in_bucket(
                key + move_pos, keysize - move_pos);
        }

        std::pair<std::string, T> getSubStr(std::pair<std::string, T>& element,
                                            size_t sub_pos = 1) {
            return std::pair<std::string, T>(element.first.substr(sub_pos),
                                             element.second);
        }

        // this func is to turn this(a hashnode) to n childs of trie_node
        // linking their child(hashnode or trienode)
        hash_node* burst(std::map<std::string, T>& elements, const CharT* key,
                         size_t keysize, trie_node* p, size_t& move_pos,
                         htrie_map* hm) {
            debugStream << "bursting with parent: " << (void*)p << std::endl;
            debugStream << "key now at " << (void*)key << " which is " << *key
                        << std::endl;
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
                    debugStream << "create new t_root: NOW THE T_ROOT IS  "
                                << (void*)hm->t_root << std::endl;
                    cur_trie_node->anode::setParent((trie_node*)hm->t_root);
                    ((trie_node*)hm->t_root)->addChildTrieNode(cur_trie_node);
                    debugStream << "root add new child: "
                                << cur_trie_node->anode::myChar << std::endl;
                    p = (trie_node*)hm->t_root;
                } else {
                    // bursting in a normal hashnode
                    cur_trie_node->anode::setParent(p);
                    p->addChildTrieNode(cur_trie_node);
                    debugStream
                        << "not-root: " << (void*)p
                        << " add new child: " << cur_trie_node->anode::myChar
                        << std::endl;
                }

                debugStream << "set new trie_node:" << (void*)cur_trie_node
                            << " with char: " << it->first << std::endl;

                if (it->first == *(key)) {
                    move_pos++;
                }

                std::map<std::string, T>& curKV = it->second;

                if (splitElements.size() == 1) {
                    return burst(curKV, key + 1, keysize - 1, cur_trie_node,
                                 move_pos, hm);
                }

                hash_node* hnode = new hash_node(burst_threshold);
                hnode->anode::setParent(cur_trie_node);
                cur_trie_node->setOnlyHashNode(hnode);
                debugStream << "create new hashnode: " << (void*)hnode
                            << std::endl;

                for (auto itt = curKV.begin(); itt != curKV.end(); itt++) {
                    std::string temp = itt->first;

                    if (temp.size() == 0) {
                        hnode->haveValue = true;
                        hnode->onlyValue = itt->second;
                        if (*(CharT*)key == cur_trie_node->anode::myChar &&
                            myTrie::hashRelative::keyEqual(temp.data(),
                                                           temp.size(), key + 1,
                                                           keysize - 1)) {
                            debugStream << ">>>>>>>>>>>>>>>find the 0-size "
                                           "target_node !!!\n ";
                            target_hashNode = hnode;
                        }
                        continue;
                    }
                    debugStream << "=============for loop=======" << std::endl;
                    debugStream << "working on string " << itt->first
                                << " with value: " << itt->second << std::endl;
                    debugStream << "curKV[i] size: " << itt->first.size()
                                << std::endl;
                    size_t hashval = myTrie::hashRelative::hash<CharT>(
                        temp.data(), temp.size());
                    debugStream << "rewrite to hnode: " << temp << std::endl;

                    // write the value to the entry
                    hnode->kvs[hashval].access_kv_in_bucket(
                        temp.data(), temp.size()) = itt->second;
                    debugStream << "move_pos is " << move_pos << std::endl;

                    if (*(CharT*)key == cur_trie_node->anode::myChar &&
                        myTrie::hashRelative::keyEqual(temp.data(), temp.size(),
                                                       key + 1, keysize - 1)) {
                        debugStream
                            << "father char compared: " << *(CharT*)key << " - "
                            << cur_trie_node->anode::myChar
                            << ">>>>>>>>>>>>>>>find the target_node !!!\n ";
                        target_hashNode = hnode;
                    }
                }

                if (hnode->need_burst()) {
                    // move 1 char
                    hash_node* thnode =
                        hnode->burst(curKV, key + 1, keysize - 1, cur_trie_node,
                                     move_pos, hm);
                    if (thnode != nullptr) {
                        target_hashNode = thnode;
                    }
                }
            }
            return target_hashNode;
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
            debugStream << "init array_bucket: " << (void*)arr_buffer << "\n";
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
            debugStream << "find in bucket: " << (void*)arr_buffer << "\n";
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
                    std::cerr << "realloc failed!!!!!!!!!1" << std::endl;
                    exit(-1);
                } else {
                    debugStream
                        << "realloc success: new_buffer:" << (void*)new_buffer
                        << std::endl;
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
            debugStream << "\n\n-----------------------------\ncall "
                           "read_array_bucket\n";
            CharT* buf = arr_buffer;
            while (!is_end_of_bucket(buf)) {
                size_t length = read_key_size(buf);
                debugStream << *((KeySizeT*)buf) << "|";

                for (size_t i = 0; i != length; i++) {
                    debugStream << *(buf + sizeof(KeySizeT) + i) << "|";
                }
                // move ptr to next header, skip keysize, string, value
                buf = buf + sizeof(KeySizeT) + length;
                debugStream << (T)*buf << "|" << std::endl;
                buf = buf + +sizeof(T);
            }
            debugStream << "END read: |" << (int)*(buf) << "| finish "
                        << std::endl;
            debugStream << "-----------------------------\n\n";
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

   public:
    anode* t_root;
    htrie_map(size_t customized_burst_threshold = DEFAULT_BURST_THRESHOLD) {
        t_root = new hash_node(customized_burst_threshold);
        burst_threshold = customized_burst_threshold;
    }

    void setRoot(anode* node) { t_root = node; }

    T& operator[](std::string key) {
        return access_operator(key.data(), key.size());
    }

    T& operator[](const CharT* key) {
        return access_operator(key, std::strlen(key));
    }

    T& access_operator(const CharT* key, size_t key_size) {
        debugStream << "-------------->accessing ";
        for (size_t i = 0; i != key_size; i++) {
            debugStream << *(key + i);
        }
        debugStream << std::endl;
        return find(key, key_size);
    }

    T& find(const CharT* key, size_t key_size) {
        anode* current_node = t_root;
        debugStream << "start finding: root is hashnode:"
                    << (t_root->isHashNode() ? "true" : "false") << " with "
                    << (void*)t_root << std::endl;

        if (t_root->isHashNode()) {
            // find a key in hash_node
            // if find: return the v reference
            // if notfind: write the key and return the v reference
            return ((hash_node*)current_node)
                ->access_kv_in_hashnode(key, key_size, this);
        }

        for (size_t pos = 0; pos < key_size; pos++) {
            if (current_node->isTrieNode()) {
                anode* parent = current_node;
                current_node =
                    ((trie_node*)current_node)->findChildNode(key[pos]);

                if (current_node == nullptr) {
                    // can't find, create a relative trie_node and a hashnode
                    // child, set current_node as the trie_node
                    current_node =
                        ((trie_node*)parent)->create_trienode_With(key[pos]);
                } else if (current_node->isHashNode()) {
                    // if find the target hashnode instead of moving the pos,
                    // pos should be recover
                    pos--;
                }
            } else {
                return ((hash_node*)current_node)
                    ->access_kv_in_hashnode(key + pos, key_size - pos, this);
            }
        }
        return ((trie_node*)current_node)
            ->getOnlyHashNode()
            ->access_onlyValue_in_hashnode(this);
    }
};  // namespace myTrie


//------------------------DEBUG, CORRECTNESS CHECK------------------------------
namespace debuging {
using std::cout;
using std::endl;
template <class CharT, class T>
void print_tree_struct(typename myTrie::htrie_map<CharT, T>::anode* root) {
    debugStream << "in node: " << (void*)root
                << " with type: " << (root->isHashNode() ? "hash" : "trie")
                << std::endl;
    // print bucket
    if (root->isHashNode()) {
        debugStream << "searching in hashnode" << endl;
        std::map<size_t, std::vector<std::pair<std::string, T>>> buckets_info =
            ((typename myTrie::htrie_map<CharT, T>::hash_node*)root)
                ->get_hash_nodes_buckets_info();
        for (auto it = buckets_info.begin(); it != buckets_info.end(); it++) {
            debugStream << it->first << ":";
            std::vector<std::pair<std::string, T>> bucket_string = it->second;
            for (int i = 0; i != bucket_string.size(); i++) {
                debugStream << bucket_string[i].first << ",";
            }
            debugStream << "|" << endl;
        }
    } else if (root->isTrieNode()) {
        std::map<CharT, typename myTrie::htrie_map<CharT, T>::anode*> childs =
            ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                ->getChildsMap();
        for (auto it = childs.begin(); it != childs.end(); it++) {
            debugStream << it->first << "\t\t";
        }
        debugStream << "\nnext level:\n";
        if (childs.size() == 0) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_struct<CharT, T>(it->second);
            }
        } else {
            debugStream << "size of child is 0" << endl;
            debugStream << "fake_root is " << (void*)root << endl;
            // print_tree_struct<CharT, T>(
            //     ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
            //         ->getOnlyHashNode());
        }
    } else {
        debugStream << "node is not trie nor hash node\n";
        exit(0);
    }
}

uint32_t trie_node_num;
uint32_t hash_node_num;
uint64_t bucket_buf_size;

template <class CharT, class T>
void print_tree_use_mem(typename myTrie::htrie_map<CharT, T>::anode* root) {
    // print bucket
    if (root->isHashNode()) {
        debugStream << "searching in hashnode" << endl;
        std::map<size_t, std::vector<std::pair<std::string, T>>> buckets_info =
            ((typename myTrie::htrie_map<CharT, T>::hash_node*)root)
                ->get_hash_nodes_buckets_info();
        for (auto it = buckets_info.begin(); it != buckets_info.end(); it++) {
            debugStream << it->first << ":";
            std::vector<std::pair<std::string, T>> bucket_string = it->second;
            for (int i = 0; i != bucket_string.size(); i++) {
                debugStream << bucket_string[i].first << ",";
            }
            debugStream << "|" << endl;
        }
    } else if (root->isTrieNode()) {
        std::map<CharT, typename myTrie::htrie_map<CharT, T>::anode*> childs =
            ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                ->getChildsMap();
        for (auto it = childs.begin(); it != childs.end(); it++) {
            debugStream << it->first << "\t\t";
        }
        debugStream << "\nnext level:\n";
        if (childs.size() == 0) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_struct<CharT, T>(it->second);
            }
        } else {
            debugStream << "size of child is 0" << endl;
            debugStream << "fake_root is " << (void*)root << endl;
            // print_tree_struct<CharT, T>(
            //     ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
            //         ->getOnlyHashNode());
        }
    } else {
        debugStream << "node is not trie nor hash node\n";
        exit(0);
    }
}

template <class CharT, class T>
void print_htrie_map(htrie_map<CharT, T> hm,
                     std::vector<std::string> checklist) {
    // correctness check
    debugStream << "finding func check:\n";
    for (auto it = checklist.begin(); it != checklist.end(); it++) {
        debugStream << "search " << *it << " got " << hm[*it] << endl;
    }

    debugStream << "-------------------------\n";

    // node check
    // print_tree_struct<CharT, T>(hm.t_root);
}

template <class CharT>
void printDiff(const CharT* key_lhs, std::size_t key_size_lhs,
               const CharT* key_rhs, std::size_t key_size_rhs, bool equal) {
    debugStream << "==comparing: \n";
    for (size_t i = 0; i != key_size_lhs; i++) {
        // debugStream << (unsigned int)*(key_lhs + i) << ',';
        debugStream << *(key_lhs + i) << ',';
    }
    debugStream << "with size:" << key_size_lhs << "\n and \n";
    for (size_t i = 0; i != key_size_rhs; i++) {
        // debugStream << (unsigned int)*(key_rhs + i) << ',';
        debugStream << *(key_rhs + i) << ',';
    }
    debugStream << "with size:" << key_size_rhs << std::endl;
    debugStream << " keyEqual res: "
                << (std::memcmp(key_lhs, key_rhs,
                                key_size_lhs * sizeof(CharT)) == 0)
                << " keyEuqal res: " << equal << std::endl;
}
}  // namespace debuging
}  // namespace myTrie

#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>
static uint64_t get_usec() {
    struct timespec tp;
    /* POSIX.1-2008: Applications should use the clock_gettime() function
       instead of the obsolescent gettimeofday() function. */
    /* NOTE: The clock_gettime() function is only available on Linux.
       The mach_absolute_time() function is an alternative on OSX. */
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return ((tp.tv_sec * 1000 * 1000) + (tp.tv_nsec / 1000));
}

#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

class NullBuffer : public std::streambuf {
   public:
    int overflow(int c) { return c; }
};

using namespace std;
int main() {
    NullBuffer null_buffer;

    streambuf* coutBuf = cout.rdbuf();
    ofstream of("out.txt");
    ofstream cc("checkAP9.txt");

    streambuf* fileBufcc = cc.rdbuf();

    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(&null_buffer);
    std::fstream f("str_normal");
    std::string url;
    uint32_t v;
    uint32_t count;
    uint64_t sta = get_usec();
    myTrie::htrie_map<char, uint32_t> hm;
    std::map<std::string, uint32_t> storeSomeStr;
    while (f >> url >> v) {
        cout << "---------------------------------------> working on num."
             << count << " url: " << url << " set value: " << v << endl;
        hm[url] = v;
        count++;
        if (count % 10000 == 0) {
            storeSomeStr[url] = v;
        }
    }
    of.flush();
    of.close();
    cout.rdbuf(coutBuf);
    f.close();
    uint64_t end = get_usec();
    cout << "finish trie_map constructing\n";

    cout << "myTrie use time: usec: " << end - sta << std::endl;

    std::fstream f4("str_normal");
    std::map<std::string, uint32_t> m;
    sta = get_usec();
    while (f4 >> url >> v) {
        m[url] = v;
    }
    end = get_usec();
    cout << "std::map use time: usec: " << end - sta << std::endl;
    f4.close();
    cout.rdbuf(&null_buffer);

    m.begin()->second = 123456;

    // checking:
    std::fstream f2("test_res_good", std::ios::out);
    std::fstream f3("test_res_wrong", std::ios::out);
    for (auto it = m.begin(); it != m.end(); it++) {
        uint32_t gotfromhm = hm[it->first];
        if (gotfromhm == it->second) {
            f2 << "good: " << it->first << " hm: " << gotfromhm
               << " map: " << it->second << std::endl;
        } else {
            f3 << "wrong answer: " << it->first << " got " << hm[it->first]
               << " from hm , actual value: " << it->second << std::endl;
            uint32_t v = it->second;
            f3 << "got from hm: ";
            for (size_t i = 0; i != sizeof(v); i++) {
                f3 << (unsigned int)*(
                          (char*)(&(hm[it->first]) + sizeof(char) * i))
                   << ",";
            }
            f3 << "\ngot from file: ";
            for (size_t i = 0; i != sizeof(v); i++) {
                f3 << (unsigned int)*((char*)(&v + sizeof(char) * i)) << ",";
            }
            f3 << std::endl;
            f3.flush();
        }
    }
    f2.close();
    f3.close();
    cout.rdbuf(coutBuf);

    std::ifstream fofwrong("test_res_wrong", std::ios::in);
    char line[1024] = {0};
    while (fofwrong.getline(line, sizeof(line))) {
        std::cerr << line << std::endl;
    }
    cout << "finish checking and printed correct/wrong res\n";

    while (true) {
        cout << "check url: ";
        cout.rdbuf(fileBufcc);
        cin >> url;
        sta = get_usec();
        cout << "check url: " << url << " got " << hm[url] << endl;
        end = get_usec();
        cout << "use " << (end - sta) / 1000 << "ms" << endl;
        cout.rdbuf(coutBuf);
    }
}