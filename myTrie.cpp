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
    // TODO: is it neccessary?
    if (key_size_lhs == 0 && key_size_rhs == 0) {
        return true;
    }
    if (key_size_lhs != key_size_rhs) {
        std::cout << "size not equal" << std::endl;
        return false;
    } else {
        return std::memcmp(key_lhs, key_rhs, key_size_lhs * sizeof(CharT)) == 0;
    }
}
template <class CharT>
void printDiff(const CharT* key_lhs, std::size_t key_size_lhs,
               const CharT* key_rhs, std::size_t key_size_rhs, bool equal) {
    std::cout << "==comparing: \n";
    for (size_t i = 0; i != key_size_lhs; i++) {
        // std::cout << (unsigned int)*(key_lhs + i) << ',';
        std::cout << *(key_lhs + i) << ',';
    }
    std::cout << "with size:" << key_size_lhs << "\n and \n";
    for (size_t i = 0; i != key_size_rhs; i++) {
        // std::cout << (unsigned int)*(key_rhs + i) << ',';
        std::cout << *(key_rhs + i) << ',';
    }
    std::cout << "with size:" << key_size_rhs << std::endl;
    std::cout << " keyEqual res: "
              << (std::memcmp(key_lhs, key_rhs, key_size_lhs * sizeof(CharT)) ==
                  0)
              << " keyEuqal res: " << equal << std::endl;
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

        hash_node* createAHashNodeWith(CharT key) {
            trie_node* new_trie_node = new trie_node(key, this);
            // TODO: remove the hard code
            new_trie_node->onlyHashNode = new hash_node(BURST_POINT);
            new_trie_node->onlyHashNode->anode::parent = new_trie_node;
            this->addChildTrieNode(new_trie_node);
            std::cout << "father: " << (void*)this
                      << " create a new trie_node: " << (void*)new_trie_node
                      << " return the onlyHahsNode: "
                      << (void*)new_trie_node->onlyHashNode << std::endl;
            return new_trie_node->onlyHashNode;
        }

        std::map<CharT, anode*> getChildsMap() { return childs; }

        anode* findChildNode(CharT c) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                std::cout << (void*)this << " have child: " << it->first
                          << " at " << (void*)it->second << std::endl;
            }
            auto found = childs.find(c);
            if (found != childs.end()) {
                return found->second;
                std::cout << "return target\n";

            } else {
                if (childs.size() == 0) {
                    std::cout << "no match child return onlyHashNode: "
                              << (void*)onlyHashNode << std::endl;
                    return onlyHashNode;
                }
                std::cout << "return nullptr\n";
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

        hash_node(size_t customized_burst_threshold) {
            anode::_node_type = node_type::HASH_NODE;
            anode::parent = nullptr;
            kvs = std::vector<array_bucket>(BUCKET_INIT_COUNT);
            burst_threshold = customized_burst_threshold;
            std::cout << "new hash_node is init, burst_threshold set to be "
                      << burst_threshold << std::endl;
            onlyValue = 999;
        }

        bool need_burst() {
            size_t node_element_count = 0;
            for (auto it = kvs.begin(); it != kvs.end(); it++) {
                node_element_count += it->bucket_element_count;
            }
            std::cout << "node element count: " << node_element_count
                      << std::endl;
            std::cout << "burst_threshold: " << burst_threshold << std::endl;
            return node_element_count > burst_threshold;
        }

        // find the key with keysize in its kvs
        T& access_kv_in_hashnode(const CharT* key, size_t keysize,
                                 htrie_map* hm) {
            hash_node* targetNode = this;
            size_t move_pos = 0;

            // TODO: burst if hash_node have too many element
            if (need_burst()) {
                std::map<std::string, T> elements;
                // add new or existed key, if existed, it will be overwritten at
                // the loop below
                elements[std::string(key, keysize)] = T{};
                for (auto it = kvs.begin(); it != kvs.end(); it++) {
                    // TODO: use the reference to reduce the extra copy
                    std::vector<std::pair<std::string, T>> bucket_element =
                        it->get_item_in_array_bucket();
                    for (auto itt = bucket_element.begin();
                         itt != bucket_element.end(); itt++) {
                        elements[itt->first] = itt->second;
                    }
                }

                for (auto it = elements.begin(); it != elements.end(); it++) {
                    std::cout << "elements: " << it->first
                              << " with value: " << it->second << std::endl;
                }
                std::cout << "bursting in " << (void*)this << std::endl;
                size_t burstPos = 0;
                // update the targetNode
                targetNode = burst(elements, key, keysize, this->anode::parent,
                                   move_pos, hm);
                std::cout << "target_node found: " << (void*)targetNode
                          << std::endl;
                delete this;
            }

            if (keysize == 0) {
                return onlyValue;
            }

            // TODO: replace the empty onlyValue
            if (targetNode->onlyValue != 999) {
                return targetNode->onlyValue;
            }

            size_t hashval = myTrie::hashRelative::hash<CharT>(
                key + move_pos, keysize - move_pos);
            std::cout << "bucket_id: " << hashval << std::endl;
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
            std::cout << "bursting with parent: " << (void*)p << std::endl;
            std::cout << "key now at " << (void*)key << " which is " << *key
                      << std::endl;
            hash_node* target_hashNode = nullptr;

            std::map<CharT, std::map<std::string, T>> splitElements;
            for (auto it = elements.begin(); it != elements.end(); it++) {
                splitElements[(it->first)[0]][it->first.substr(1)] = it->second;
            }

            // TODO: improve
            if (splitElements.size() == 1) {
                CharT k = splitElements.begin()->first;
                std::map<std::string, T> newElement =
                    splitElements.begin()->second;
                trie_node* cur_trie_node = new trie_node(k, nullptr);

                if (p == nullptr) {
                    // bursting in a root hashnode
                    // the t_root is update to a empty trie_node
                    hm->t_root = new trie_node('\0', nullptr);
                    std::cout << "create new t_root: NOW THE T_ROOT IS  "
                              << (void*)hm->t_root << std::endl;
                    cur_trie_node->anode::setParent((trie_node*)hm->t_root);
                    ((trie_node*)hm->t_root)->addChildTrieNode(cur_trie_node);
                    std::cout << "root add new child: "
                              << cur_trie_node->anode::myChar << std::endl;
                    p = (trie_node*)hm->t_root;
                } else {
                    // bursting in a normal hashnode
                    cur_trie_node->anode::setParent(p);
                    p->addChildTrieNode(cur_trie_node);
                    std::cout << "not-root: " << (void*)p << " add new child: "
                              << cur_trie_node->anode::myChar << std::endl;
                }

                if (k == *(key)) {
                    move_pos++;
                }

                return burst(newElement, key + 1, keysize - 1, cur_trie_node,
                             move_pos, hm);
            }

            for (auto it = splitElements.begin(); it != splitElements.end();
                 it++) {
                trie_node* cur_trie_node = new trie_node(it->first, nullptr);

                if (p == nullptr) {
                    // bursting in a root hashnode
                    // the t_root is update to a empty trie_node
                    hm->t_root = new trie_node('\0', nullptr);
                    std::cout << "create new t_root: NOW THE T_ROOT IS  "
                              << (void*)hm->t_root << std::endl;
                    cur_trie_node->anode::setParent((trie_node*)hm->t_root);
                    ((trie_node*)hm->t_root)->addChildTrieNode(cur_trie_node);
                    std::cout << "root add new child: "
                              << cur_trie_node->anode::myChar << std::endl;
                    p = (trie_node*)hm->t_root;
                } else {
                    // bursting in a normal hashnode
                    cur_trie_node->anode::setParent(p);
                    p->addChildTrieNode(cur_trie_node);
                    std::cout << "not-root: " << (void*)p << " add new child: "
                              << cur_trie_node->anode::myChar << std::endl;
                }

                std::cout << "set new trie_node:" << (void*)cur_trie_node
                          << " with char: " << it->first << std::endl;

                // TODO: move_pos
                if (it->first == *(key)) {
                    move_pos++;
                }

                std::map<std::string, T>& curKV = it->second;
                hash_node* hnode = new hash_node(burst_threshold);
                hnode->anode::setParent(cur_trie_node);
                cur_trie_node->setOnlyHashNode(hnode);
                std::cout << "create new hashnode: " << (void*)hnode
                          << std::endl;

                for (auto itt = curKV.begin(); itt != curKV.end(); itt++) {
                    std::string temp = itt->first;

                    if (temp.size() == 0) {
                        hnode->onlyValue = itt->second;
                        if (*(CharT*)key == cur_trie_node->anode::myChar &&
                            myTrie::hashRelative::keyEqual(temp.data(),
                                                           temp.size(), key + 1,
                                                           keysize - 1)) {
                            std::cout << ">>>>>>>>>>>>>>>find the 0-size "
                                         "target_node !!!\n ";
                            target_hashNode = hnode;
                        }
                        continue;
                    }
                    std::cout << "=============for loop=======" << std::endl;
                    std::cout << "working on string " << itt->first
                              << " with value: " << itt->second << std::endl;
                    std::cout << "curKV[i] size: " << itt->first.size()
                              << std::endl;
                    size_t hashval = myTrie::hashRelative::hash<CharT>(
                        temp.data(), temp.size());
                    std::cout << "rewrite to hnode: " << temp << std::endl;

                    // write the value to the entry
                    hnode->kvs[hashval].access_kv_in_bucket(
                        temp.data(), temp.size()) = itt->second;
                    std::cout << "move_pos is " << move_pos << std::endl;

                    // if (myTrie::hashRelative::keyEqual(temp.data(),
                    // temp.size(),
                    //                                    key, keysize)) {
                    //     myTrie::hashRelative::printDiff(
                    //         temp.data(), temp.size(), key, keysize,
                    //         true);

                    //     std::cout << ">>>>>>>>>>>>>>>find the
                    //     target_node!!!\n"; target_hashNode = hnode;
                    // } else {
                    //     myTrie::hashRelative::printDiff(
                    //         temp.data(), temp.size(), key, keysize,
                    //         false);
                    // }

                    if (*(CharT*)key == cur_trie_node->anode::myChar &&
                        myTrie::hashRelative::keyEqual(temp.data(), temp.size(),
                                                       key + 1, keysize - 1)) {
                        // myTrie::hashRelative::printDiff(temp.data(),
                        //                                 temp.size(), key
                        //                                 + 1, keysize - 1,
                        //                                 true);

                        std::cout
                            << "father char compared: " << *(CharT*)key << " - "
                            << cur_trie_node->anode::myChar
                            << ">>>>>>>>>>>>>>>find the target_node !!!\n ";
                        target_hashNode = hnode;
                    } else {
                        // myTrie::hashRelative::printDiff(temp.data(),
                        //                                 temp.size(), key
                        //                                 + 1, keysize - 1,
                        //                                 false);
                    }

                    // if (myTrie::hashRelative::keyEqual(temp.data(),
                    // temp.size(),
                    //                                    key + move_pos,
                    //                                    keysize -
                    //                                    move_pos)) {
                    //     myTrie::hashRelative::printDiff(
                    //         temp.data(), temp.size(), key + move_pos,
                    //         keysize - move_pos, true);

                    //     std::cout
                    //         << ">>>>>>>>>>>>>>>find the target_node !!!\n
                    //         ";
                    //     target_hashNode = hnode;
                    // } else {
                    //     myTrie::hashRelative::printDiff(
                    //         temp.data(), temp.size(), key + move_pos,
                    //         keysize - move_pos, false);
                    // }
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
            std::cout << "init array_bucket: " << (void*)arr_buffer << "\n";
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
            std::cout << "find in bucket: " << (void*)arr_buffer << "\n";
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
                    std::cout << "realloc failed!!!!!!!!!1" << std::endl;
                    exit(-1);
                } else {
                    std::cout
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
        std::cout << "-------------->accessing ";
        for (size_t i = 0; i != key_size; i++) {
            std::cout << *(key + i);
        }
        std::cout << std::endl;
        return find(key, key_size);
    }

    T& find(const CharT* key, size_t key_size) {
        anode* current_node = t_root;
        std::cout << "start finding: root is hashnode:"
                  << (t_root->isHashNode() ? "true" : "false") << " with "
                  << (void*)t_root << std::endl;
        for (size_t pos = 0; pos < key_size; pos++) {
            // didn't reach the leaf
            if (current_node->isTrieNode()) {
                std::cout << "current_node is trienode: " << (void*)current_node
                          << " with char: " << current_node->anode::myChar
                          << std::endl;
                anode* father_node = current_node;
                current_node =
                    ((trie_node*)current_node)->findChildNode(key[pos]);

                // can't find, need to create a trie_node with a hashtable
                // child
                if (current_node == nullptr) {
                    current_node = ((trie_node*)father_node)
                                       ->createAHashNodeWith(key[pos]);
                    pos++;
                    std::cout
                        << "--current_node is hashnode: " << (void*)current_node
                        << std::endl;
                    std::cout << "pos:" << pos
                              << " key_size - pos: " << key_size - pos
                              << " looking for ";
                    for (size_t i = 0; i != key_size - pos; i++) {
                        std::cout << *(key + i + pos);
                    }
                    std::cout << std::endl;

                    // find a key in hash_node
                    // if find: return the v reference
                    // if notfind: write the key and return the v reference
                    return ((hash_node*)current_node)
                        ->access_kv_in_hashnode(key + pos, key_size - pos,
                                                this);
                } else {
                    std::cout << "find chilld res: " << (void*)current_node
                              << " is a "
                              << (current_node->isHashNode() ? "hash" : "trie")
                              << std::endl;
                    std::cout << "pos: " << pos << " keysize: " << key_size
                              << std::endl;
                    if (current_node->isHashNode()) {
                        std::cout
                            << "-current_node is hashnode: dir access! only T"
                            << ((hash_node*)current_node)->onlyValue
                            << std::endl;
                        std::cout << "-current_node is hashnode: "
                                  << (void*)current_node << std::endl;
                        std::cout << "pos:" << pos
                                  << " key_size - pos: " << key_size - pos
                                  << " looking for ";
                        for (size_t i = 0; i != key_size - pos; i++) {
                            std::cout << *(key + i + pos);
                        }
                        std::cout << std::endl;
                        return ((hash_node*)current_node)
                            ->access_kv_in_hashnode(key + pos, key_size - pos,
                                                    this);
                    } else {
                    }
                }
                // else if (current_node->isHashNode()) {
                //     std::cout << "-current_node is hashnode" <<
                //     std::endl; std::cout << "looking for "; for (size_t i
                //     = 0; i != key_size; i++) {
                //         std::cout << *(key + i + pos);
                //     }
                //     std::cout << std::endl;
                //     return ((hash_node*)current_node)
                //         ->access_kv_in_hashnode(key + pos, key_size -
                //         pos,
                //                                 this);
                // }

            } else {
                std::cout << "---current_node is hashnode: "
                          << (void*)current_node << std::endl;
                std::cout << "pos:" << pos
                          << " key_size - pos: " << key_size - pos
                          << " looking for ";
                for (size_t i = 0; i != key_size - pos; i++) {
                    std::cout << *(key + i + pos);
                }
                std::cout << std::endl;

                if (pos == 0) {
                    // find a key in hash_node
                    // if find: return the v reference
                    // if notfind: write the key and return the v reference
                    return ((hash_node*)current_node)
                        ->access_kv_in_hashnode(key + pos, key_size - pos,
                                                this);
                } else {
                    return ((hash_node*)current_node)
                        ->access_kv_in_hashnode(key + pos - 1,
                                                key_size - pos + 1, this);
                }
            }
        }
        std::cout << "return the only value" << std::endl;
        return ((trie_node*)current_node)->getOnlyHashNode()->onlyValue;
    }
};

using std::cout;
using std::endl;
template <class CharT, class T>
void print_tree_struct(typename myTrie::htrie_map<CharT, T>::anode* root) {
    cout << "in node: " << (void*)root
         << " with type: " << (root->isHashNode() ? "hash" : "trie")
         << std::endl;
    // print bucket
    if (root->isHashNode()) {
        std::cout << "searching in hashnode" << endl;
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
        if (childs.size() == 0) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_struct<CharT, T>(it->second);
            }
        } else {
            cout << "size of child is 0" << endl;
            cout << "fake_root is " << (void*)root << endl;
            // print_tree_struct<CharT, T>(
            //     ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
            //         ->getOnlyHashNode());
        }
    } else {
        cout << "node is not trie nor hash node\n";
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
        std::cout << "searching in hashnode" << endl;
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
        if (childs.size() == 0) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_struct<CharT, T>(it->second);
            }
        } else {
            cout << "size of child is 0" << endl;
            cout << "fake_root is " << (void*)root << endl;
            // print_tree_struct<CharT, T>(
            //     ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
            //         ->getOnlyHashNode());
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
    // print_tree_struct<CharT, T>(hm.t_root);
}
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
using namespace std;
int main() {
    // using namespace std;
    // using namespace myTrie;

    // vector<string> tests;
    // tests.push_back("abc");
    // tests.push_back("abcd");
    // tests.push_back("bbcd");
    // tests.push_back("bcde");

    // tests.push_back("bcdef");
    // tests.push_back("bcded");
    // tests.push_back("bcabc");
    // tests.push_back("ccc");

    // htrie_map<char, int> hm;
    // for (auto it = tests.begin(); it != tests.end(); it++)
    //     hm[*it] = it - tests.begin();
    // std::cout << "------------------------\n";
    // // print_htrie_map<char, int>(hm, tests);
    // while (true) {
    //     std::string w;
    //     cin >> w;
    //     cout << "search " << w << " got " << hm[w] << endl;
    // }
    streambuf* coutBuf = cout.rdbuf();
    ofstream of("out.txt");
    ofstream cc("checkAP9.txt");

    streambuf* fileBufcc = cc.rdbuf();

    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
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

    std::cout << "myTrie use time: usec: " << end - sta << std::endl;
    std::fstream f4("str_normal");
    std::map<std::string, uint32_t> m;
    sta = get_usec();
    while (f4 >> url >> v) {
        m[url] = v;
    }
    end = get_usec();
    std::cout << "std::map use time: usec: " << end - sta << std::endl;
    f4.close();

    // checking:
    std::fstream f2("test_res_good", std::ios::out);
    std::fstream f3("test_res_wrong", std::ios::out);
    for (auto it = m.begin(); it != m.end(); it++) {
        if (hm[it->first] == it->second) {
            f2 << "good: " << it->first << std::endl;
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

    std::cout << "finish checking and printed correct/wrong res\n";

    while (true) {
        cout.rdbuf(fileBufcc);
        cout << "check url: ";
        cin >> url;
        sta = get_usec();
        cout << "check url: " << url << " got " << hm[url] << endl;
        end = get_usec();
        cout << "use " << (end - sta) / 1000 << "ms" << endl;
        cout.rdbuf(coutBuf);
    }
}