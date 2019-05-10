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
#include "../util/my_timer.hpp"

#include <fstream>
#include <ios>
#include <iostream>
#include <set>
#include <string>

#include <assert.h>

//for min(x,y)
#include <algorithm>

#include <boost/unordered_map.hpp>

/* 
 * DEFAULT_Associativity * DEFAULT_Bucket_num should be set to be greater
 * than 26 for 26 alaphbet and ", @ .etc.(the test in lubm40 shows that it
 * has 50 char species)
 */
#define DEFAULT_Associativity 8
#define DEFAULT_Bucket_num 10
#define DEFAULT_Max_bytes_per_kv 4096
#define DEFAULT_SPECIAL_Max_bytes_per_kv_RATIO 4
#define DEFAULT_Burst_ratio 0.75

#define ALPHABET 256

#define DEFAULT_SPECIAL_Max_bytes_per_kv \
    DEFAULT_Max_bytes_per_kv* DEFAULT_SPECIAL_Max_bytes_per_kv_RATIO

#define NBITS_SPECIAL 1
#define NBITS_LEN 7   // 128 length
#define NBITS_POS 12  // 4096(1 align) pos
#define NBITS_PID 12  // 4096 page

#define NBITS_SPECIAL_S 1
#define NBITS_LEN_S 13  // 8192 length
#define NBITS_POS_S 9   // 512(32 align) pos
#define NBITS_PID_S 9   // 512 page

#define MAX_NORMAL_PAGE (1 << NBITS_PID)
#define MAX_SPECIAL_PAGE (1 << NBITS_PID_S)

#define MAX_NORMAL_LEN (1 << NBITS_LEN)

// #define DEFAULT_SPECIAL_ALIGNMENT \
//     DEFAULT_SPECIAL_Max_bytes_per_kv / (1 << NBITS_POS_S)
#define DEFAULT_SPECIAL_ALIGNMENT 32
#define DEFAULT_NORMAL_ALIGNMENT 1

#define FAST_PATH_NODE_NUM (20)

static uint32_t longest_string_size;

uint32_t recal_element_num_of_1st_char_counter = 0;
uint32_t burst_total_counter = 0;

uint64_t burst_total_time = 0;
uint64_t cal_prefix_total_time = 0;
uint64_t write_page_total_time = 0;

// todo: wait to be deleted, just for recording the time that expand() cost
uint64_t expand_cost_time = 0;
uint64_t rehash_cost_time = 0;
uint64_t rehash_total_num = 0;

uint64_t shrink_total_time = 0;
uint64_t clean_prefix_total_time = 0;

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
    public:
    
    /* trie node type */
    class trie_node;
    class hash_node;

    /* iterator that contains */
    class iterator;

    class slot;
    /* page management */
    class page_manager;
    class page_group;
    class page;

    /* group type divided by the keysize */
    enum class group_type : unsigned char {
        NORMAL_GROUP,
        SPECIAL_GROUP,
    };

    /* helper function */
    static inline group_type get_group_type(size_t keysize) {
        return keysize < MAX_NORMAL_LEN ? group_type::NORMAL_GROUP
                                        : group_type::SPECIAL_GROUP;
    }

    static inline group_type get_group_type(slot* s) {
        return s->get_length() < MAX_NORMAL_LEN ? group_type::NORMAL_GROUP
                                                : group_type::SPECIAL_GROUP;
    }

    /* helper class for hash_node get its page_groups */
    class page_group_package {
        typename page_manager::page_group* n_group;
        typename page_manager::page_group* s_group;

       public:
        page_group_package(typename page_manager::page_group* ng,
                           typename page_manager::page_group* sg)
            : n_group(ng), s_group(sg) {}

        inline typename page_manager::page_group* get_page_group(slot* s) {
            return get_group_type(s->get_length()) == group_type::SPECIAL_GROUP
                       ? s_group
                       : n_group;
        }

        inline typename page_manager::page_group* get_page_group(group_type get_type) {
            return get_type == group_type::SPECIAL_GROUP ? s_group : n_group;
        }

        inline void set_page_group(
            group_type get_type,
                typename page_manager::page_group* update_page_group) {
          get_type == group_type::SPECIAL_GROUP ? (s_group = update_page_group)
                                                : (n_group = update_page_group);
        }

        inline bool is_special(slot* s) {
            return get_group_type(s->get_length()) == group_type::SPECIAL_GROUP;
        }

        // get function
        inline char* get_content_pointer(slot* s) {
            return is_special(s) ? s_group->get_content_pointer_in_page(s)
                                 : n_group->get_content_pointer_in_page(s);
        }

        inline T get_value(slot* s) {
            return is_special(s) ? s_group->get_value_in_page(s)
                                 : n_group->get_value_in_page(s);
        }
    };

    /* node type definition */
    enum class node_type : unsigned char { HASH_NODE, TRIE_NODE };
    /* node's base class */
    class anode {
       private:
        node_type n_type_;
        trie_node* parent_;

        bool have_value_;
        T value_;

       public:
       // TODO: adopt a more save-memory way 
        string prefix_;

       public:
        anode(node_type n_type, trie_node* parent, const CharT *key, size_t key_size)
            : n_type_(n_type),
              parent_(parent),
              have_value_(false),
              value_(T()),
              prefix_("") {
            // If current node's layer equals to multiple of FAST_PATH_NODE_NUM,
            // we set up a fast path in its x-nd grandparent
            if (key_size % FAST_PATH_NODE_NUM == 0 && key_size != 0)
                set_up_fast_path(key + key_size - FAST_PATH_NODE_NUM,
                                        FAST_PATH_NODE_NUM);
        }

        bool is_hash_node() { return n_type_ == node_type::HASH_NODE; }
        bool is_trie_node() { return n_type_ == node_type::TRIE_NODE; }

        void set_parent(trie_node* p) { parent_ = p; }
        trie_node* get_parent() { return parent_; }

        node_type get_node_type() { return n_type_; }

        // virtual function for page_manager resize
        virtual void traverse_for_pgm_resize(page_manager* old_pm,
                                             page_manager* new_pm,
                                             group_type resize_type) = 0;

        // for level traverse
        virtual vector<anode*> print_node_info() = 0;

        iterator insert_value_in_node(const string &prefix, T v,
                                      htrie_map<CharT, T>* hm) {
            value_ = v;
            have_value_ = true;
            hm->set_v2k(v, this, -1);
            prefix_ = prefix;
            return iterator(have_value_, value_, this, -1, -1);
        }

        iterator search_kv_in_node() {
            return iterator(have_value_, value_, this, -1, -1);
        }

        void set_prefix(const string& prefix) { prefix_ = prefix; }
        const string& get_prefix() { return prefix_; }

        // Get the x-nd grand parent for adding a fast path
        inline trie_node* get_fast_path_parent() {
            trie_node* cur_parent = (trie_node*)this;
            for (int i = 0; i != FAST_PATH_NODE_NUM; i++) {
                cur_parent = cur_parent->anode::get_parent();
                if (cur_parent == nullptr) return nullptr;
            }
            return cur_parent;
        }

        // Set up a fast path in its target grandparent
        void set_up_fast_path(const char* key, size_t key_size) {
            // Get the target parent who is going to add a fast path for destination node(this)
            trie_node* add_fast_path_parent = get_fast_path_parent();
            assert(add_fast_path_parent != nullptr);
            add_fast_path_parent->add_fast_path(key, key_size, this);
        }
    };

    /*-------------fast_path_manager-------------*/
    class fast_path_manager {
       private:
        struct fast_path {
            unsigned int hash_val_;
            // TODO: store in page_manager
            string fast_path_string_;
            anode* dest_node_;

           public:
            fast_path(unsigned int hash_val, const string &fast_path_string, anode* dest_node)
                : hash_val_(hash_val),
                  fast_path_string_(fast_path_string),
                  dest_node_(dest_node) {}

            inline void set_dest_node(anode *node) { dest_node_ = node; }

            inline unsigned int get_hash_val() { return hash_val_; }
            inline string& get_string() { return fast_path_string_; }
            inline anode* get_dest_node() { return dest_node_; }
        };

        vector<fast_path> fast_paths;

       public:
        fast_path_manager() {}

        inline void insert_fast_path(const char* key, size_t key_size, anode* node_ptr) {
            //Insert new element
            fast_path new_fast_path(
                myTrie::hashRelative::hash(key, key_size, 1), string(key, key_size), node_ptr);

            for (auto it = fast_paths.begin(); it != fast_paths.end();
                 it++) {
              if (it->get_hash_val() > new_fast_path.get_hash_val()) {
                fast_paths.insert(it, new_fast_path);
                return;
              } else if (it->get_hash_val() == new_fast_path.get_hash_val()) {
                // TODO: FIXME: the hash value may be the same
                // Maybe the list to solve the same hash_val, different string problem(collision)
                it->set_dest_node(node_ptr);
                return;
              }
            }
            fast_paths.push_back(new_fast_path);
        }

        // TODO: FIXME: consider the same hash value situation: consider the list
        // | 5 | 5 | 5 | 5 | 5 | 5 | 5 |
        inline anode* lookup_fast_path(const char* key,
                                                 size_t key_size) {
            // binary search
            unsigned int target_hash_val =
                myTrie::hashRelative::hash(key, key_size);

            size_t node_size = fast_paths.size();
            int low = 0;
            int high = node_size - 1;
            while (low < high) {
                int mid = (low + high) >> 1;
                if (fast_paths[mid].get_hash_val() < target_hash_val) {
                    low = mid + 1;
                } else
                    high = mid;
            }

            // check the same hash value situation
            for (int i = low;
                 low != node_size &&
                 fast_paths[i].get_hash_val() == target_hash_val;
                 i++) {
                if (myTrie::hashRelative::keyEqual(
                        fast_paths[i].get_string().data(),
                        fast_paths[i].get_string().size(), key,
                        key_size)) {
                    return fast_paths[i].get_dest_node();
                }
            }
            return nullptr;
        }

        inline size_t size() { return fast_paths.size(); }

        unsigned int get_fpm_memory() {
            return size() * (sizeof(fast_path) + FAST_PATH_NODE_NUM);
        }
    };

#ifdef ARRAY_REP
    /*-----------------array child_representation-----------------*/
    class child_representation {
        trie_node** childs;
        int number;

       public:
        class child_iterator {
           public:
            char first;
            trie_node* second;

            child_iterator(char c, trie_node* tn) : first(c), second(tn) {}

            inline child_iterator* operator->() { return this; }

            inline bool operator==(child_iterator& right) {
                return second == right.second;
            }

            inline bool operator!=(const child_iterator& right) {
                return second != right.second;
            }
        };

        child_representation() : number(0) {
            childs = (trie_node**)malloc(sizeof(trie_node*) * ALPHABET);
            for (int i = 0; i != ALPHABET; i++) {
                childs[i] = nullptr;
            }
        }

        inline trie_node*& operator[](char c) {
            number++;
            return childs[(int)c];
        }

        inline child_iterator find(char c) {
            return child_iterator(c, childs[(int)c]);
        }

        inline size_t size() { return number; }

        inline child_iterator end() { return child_iterator(0, nullptr); }

        inline child_iterator begin() {
            if (number != 0) {
                for (int i = 0; i != ALPHABET; i++) {
                    if (childs[i]) return child_iterator((char)i, childs[i]);
                }
            }
            return child_iterator(0, nullptr);
        }

        inline void get_childs(vector<anode*>& res) {
            for (int i = 0; i != ALPHABET; i++) {
                if (childs[i]) res.push_back(childs[i]);
            }
        }

        inline void get_childs_with_char(vector<pair<CharT, trie_node*>>& res) {
            for (int i = 0; i != ALPHABET; i++)
                if (childs[i])
                    res.push_back(pair<CharT, trie_node*>((char)i, childs[i]));
        }

        size_t get_childs_representation_mem() {
            return sizeof(child_representation) + ALPHABET * sizeof(trie_node*);
        }
    };
#endif

#ifdef LIST_REP
    /*-----------------list child_representation-----------------*/
    class child_representation {
        class child_node {
           public:
            char child_node_char;
            anode* current;
            child_node* next;

            child_node(char cnc, anode* cur)
                : child_node_char(cnc), current(cur), next(nullptr) {}

            child_node()
                : child_node_char(0), current(nullptr), next(nullptr) {}

            inline bool have_next() { return next != nullptr; }

            inline child_node* next_child() { return next; }

            inline void add_next_child(char c) {
                next = new child_node(c, nullptr);
            }
        };

        child_node* first_child_;
        int size_;

       public:
        child_representation() : size_(0), first_child_(nullptr) {}

        inline anode*& operator[](char c) {
            //find the ok node
            child_node* current_child_node = first_child_;
            child_node* last_child_node = nullptr;

            while (current_child_node) {
                if (current_child_node->child_node_char == c)
                    return current_child_node->current;

                last_child_node = current_child_node;
                current_child_node = current_child_node->next;
            }

            // find no target node, add one
            if (first_child_ == nullptr) {
                first_child_ = new child_node(c, nullptr);

                size_++;
                return first_child_->current;
            }

            last_child_node->add_next_child(c);
            size_++;
            return last_child_node->next_child()->current;
        }

        inline anode* find(char c) {
            child_node* current_child_node = first_child_;
            while(current_child_node != nullptr) {
                if (current_child_node->child_node_char == c)
                    return current_child_node->current;

                current_child_node = current_child_node->next_child();
            };
            return nullptr;
        }

        inline size_t size() { return size_; }

        /* for traverse_for_pgm_resize() */
        inline void get_childs(vector<anode*>& res) {
            child_node* current_child_node = first_child_;

            while (current_child_node) {
                res.push_back(current_child_node->current);
                current_child_node = current_child_node->next;
            }
        }
        
        // TODO
        /* for shrink node */
        inline void get_childs_with_char(vector<pair<CharT, anode*>>& res) {
            child_node* current_child_node = first_child_;

            while (current_child_node) {
                res.push_back(
                    pair<CharT, anode*>(current_child_node->child_node_char,
                                            current_child_node->current));

                current_child_node = current_child_node->next;
            }
        }

        ~child_representation() {
            // release the list
            child_node* current_child_node = first_child_;
            child_node* previous_current_child_node = nullptr;

            while (current_child_node) {
                previous_current_child_node = current_child_node;
                current_child_node = current_child_node->next;
                delete (previous_current_child_node);
            }
        }

        /* helper function: memory evaluation */
        size_t get_childs_representation_mem() {
            return sizeof(child_representation) + size_ * sizeof(child_node);
        }
    };
#endif

#ifdef MAP_REP
    /*-----------------map child_representation-----------------*/
    class child_representation {

       public:
        map<char, trie_node*> childs;
        child_representation() {}

        inline trie_node*& operator[](char c) { return childs[(int)c]; }

        inline typename map<char, trie_node*>::iterator find(char c) {
            return childs.find(c);
        }

        inline size_t size() { return childs.size(); }

        inline typename map<char, trie_node*>::iterator end() {
            return childs.end();
        }

        inline typename map<char, trie_node*>::iterator begin() {
            return childs.begin();
        }

        void get_childs(vector<anode*> &res) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                res.push_back(it->second);
            }
        }

        inline void get_childs_with_char(vector<pair<CharT, trie_node*>>& res) {
            for(auto it = childs.begin();it!=childs.end();it++){
                res.push_back(pair<CharT,trie_node*>(it->first, it->second));
            }
        }
    };
#endif

    class trie_node : public anode {
       private:
        friend void htrie_map::traverse_level();
        fast_path_manager *fpm_;
        child_representation childs_;  // store the suffix of hash_node or trie_node

       public:
        trie_node(trie_node* p, const char* key, size_t key_size)
            : anode(node_type::TRIE_NODE, p, key, key_size), fpm_(nullptr) {}

        // add a fast path of string(key, key_size) in fast path manager
        void add_fast_path(const char* key, size_t key_size, anode* node) {
            if (fpm_ == nullptr)
                fpm_ = new fast_path_manager();
            fpm_->insert_fast_path(key, key_size, node);
            return;
        }

        // TODO: deconstructor

        // the virtual function for page_manager page resize
        void traverse_for_pgm_resize(page_manager* old_pm, page_manager* new_pm,
                                     group_type resize_type) {
            vector<anode*> childs;
            get_childs_vector(childs);
            for (auto it : childs) {
                it->traverse_for_pgm_resize(old_pm, new_pm, resize_type);
            }
        }

        vector<anode*> print_node_info(){
            vector<anode*> res;
            get_childs_vector(res);
            cout << "t:" << (void*)this << " ";

            vector<pair<CharT, anode*>> res2;
            childs_.get_childs_with_char(res2);
            for (auto c : res2)
                cout << ((c.second)->is_hash_node() ? "h:" : "t:") << c.first
                     << " " << c.second << "  ";
            return res;
        }

        /*---helper function for traverse_for_pgm_resize----*/
        void get_childs_vector(vector<anode*> &res) {
            //get all trie_node childs from child_representation
            childs_.get_childs(res);
        }

        // add node in child representation
        void add_child(const CharT c, anode* node) { childs_[c] = node; }

        // finding target, if target doesn't exist, create new trie_node with
        // hash_node son and return new trie_node
        anode* find_trie_node_child(bool findMode, const CharT* key, size_t &pos,
                                    size_t key_size, htrie_map<CharT, T>* hm) {
            // Find in fast path
            // If find the target anode in fpm(fast path manager), we return the
            // fast_path_node
            if (fpm_ != nullptr && (pos + FAST_PATH_NODE_NUM < key_size)) {
              anode *fast_path_node =
                  fpm_->lookup_fast_path(key + pos, FAST_PATH_NODE_NUM);
              if (fast_path_node != nullptr) {
                pos += FAST_PATH_NODE_NUM;
                return fast_path_node;
              }
            }

            // Find in normal path
            anode* target_node = childs_.find(key[pos]);
            pos++;
            if (findMode) {
                return target_node;
            } else {
                return target_node;
            }
        }
    };

    class hash_node : public anode {
       public:
        slot* key_metas;
        size_t elem_num;

        size_t cur_associativity = 1;

        // normal page_group id
        uint8_t normal_pgid;
        // special page_group id
        uint8_t special_pgid;

       public:
        // debug function
        void print_slot(int i, int j, htrie_map<CharT, T>* hm) {
            slot* s = key_metas + i * cur_associativity + j;
            cout << i * cur_associativity + j << ":" << s->get_special() << ","
                 << s->get_length() << "," << s->get_pos() << ","
                 << s->get_page_id() << ",";
            string str = string(hm->get_content_pointer(get_page_group_id(s), s),
                                s->get_length());
            cout << str;
            T v = hm->get_value(get_page_group_id(s), s);
            cout << "=" << v << "\n";
        }

        static void print_slot(slot* s) {
            cout << s->get_special() << "," << s->get_length() << ","
                 << s->get_pos() << "," << s->get_page_id() << "," << endl;
        }

        void print_key_metas(htrie_map<CharT, T>* hm) {
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != cur_associativity; j++) {
                    print_slot(i, j, hm);
                }
                cout << "---\n";
            }
            cout << endl;
        }

       public:
        explicit hash_node(trie_node* p,const string &prefix, htrie_map<CharT, T>* hm,
                           size_t need_associativity = 1)
            : anode(node_type::HASH_NODE, p, prefix.data(), prefix.size()),
              cur_associativity(need_associativity > Associativity
                                    ? Associativity
                                    : need_associativity),
              elem_num(0),
              normal_pgid(hm->get_normal_group_id()),
              special_pgid(hm->get_special_group_id()) {
            anode::set_prefix(prefix);

            key_metas = new slot[cur_associativity * Bucket_num]();
        }

        void traverse_for_pgm_resize(page_manager* old_pm, page_manager* new_pm,
                                     group_type resize_type) {
            size_t resize_pgid = -1;
            size_t new_pgid = -1;
            if (resize_type == group_type::SPECIAL_GROUP) {
                resize_pgid = special_pgid;
                special_pgid = new_pm->require_group_id(group_type::SPECIAL_GROUP);
                new_pgid = special_pgid;
            } else {
                resize_pgid = normal_pgid;
                normal_pgid = new_pm->require_group_id(group_type::NORMAL_GROUP);
                new_pgid = normal_pgid;
            }

            assert(resize_pgid != -1 && new_pgid != -1);

            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != cur_associativity; j++) {
                    slot* s = get_slot(i, j);

                    // ignore the slot that not belong to current resize group
                    // type
                    if (s->is_empty() || get_group_type(s) != resize_type)
                        continue;

                    // get the content from old page_manager and write it to the
                    // new page_manager
                    s->set_slot(new_pm->write_kv(
                        new_pgid,
                        old_pm->get_content_pointer_in_pm(resize_pgid, s),
                        s->get_length(),
                        old_pm->get_value_in_pm(resize_pgid, s)));
                }
            }
        }

        vector<anode*> print_node_info(){
            cout << "h:" << elem_num << " " << (void*)this << "  ";
            return vector<anode*>();
        }

        ~hash_node() { delete[] key_metas; }

        inline slot* get_slot(size_t bucketid, size_t slotid) {
            return key_metas + bucketid * cur_associativity + slotid;
        }

        inline slot* get_slot(int index) { return key_metas + index; }

        inline int get_index(slot* s) { return s - key_metas; }

        inline int get_index(int bucketid, int slotid) {
            return bucketid * cur_associativity + slotid;
        }

        /* 
         * For eliminating the index update in expand_key_metas_space
         * we store the column-store-index in v2k instead of row-store-index
         */
        inline slot* get_column_store_slot(int column_store_index) {
            return key_metas +
                   (cur_associativity * (column_store_index % Bucket_num)) +
                   column_store_index / Bucket_num;
        }

        inline int get_column_store_index(slot* s) {
            return Bucket_num * (get_index(s) % cur_associativity) +
                   (get_index(s) / cur_associativity);
        }

        /*---------function that changed key_metas layout---------*/

        /*------------------ 0. helper function------------------*/

        /*------------------ 1. expand function------------------*/
        bool expand_key_metas_space() {
            uint64_t sta = get_time();

            // Already max associativity
            // We cannot expand anymore, return false
            if (cur_associativity == Associativity) {
                return false;
            }

            // Get the associativity we need, expand 2 times of cur_associativity
            unsigned int need_associativity = cur_associativity << 1;
            if (need_associativity > Associativity) {
                need_associativity = Associativity;
            }

            // Allocate a bigger memory for new key_metas
            slot* new_key_metas = new slot[need_associativity * Bucket_num]();

            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != need_associativity; j++) {
                    slot* cur_new_slot =
                        new_key_metas + i * need_associativity + j;
                    if (j < cur_associativity) {
                        slot* cur_slot = key_metas + i * cur_associativity + j;
                        cur_new_slot->set_slot(cur_slot);
                    } else {
                        cur_new_slot->set_slot(0, 0, 0, 0);
                    }
                }
            }

            // Switch the old key_metas to the new key_metas and release the old
            // key_metas
            delete[] key_metas;
            key_metas = new_key_metas;

            // update current associativity
            cur_associativity = need_associativity;
            uint64_t end = get_time();

            expand_cost_time += end - sta;

            return true;
        }

        /*------------------ 2. cuckoo hash function------------------*/
        /*  cuckoo hash helper function  */
        // Return previous slotid in current bucket
        // If current slotid is 0, return -1
        inline int get_previous_slotid_in_same_bucket(int slotid) {
            return slotid == 0 ? -1 : slotid - 1;
        }

        // Return the first empty slot or last slot in current bucket
        inline int get_last_slotid_in_bucket(int bucketid) {
            for (int i = 0; i != cur_associativity; i++)
                if (get_slot(bucketid, i)->is_empty()) return i;

            return cur_associativity - 1;
        }

        // Return another possible bucketid that the slot *s can be at
        inline size_t get_another_bucketid(page_group_package& pgp, slot* s,
                                           size_t current_bucketid) {
            const char* key = pgp.get_content_pointer(s);
            size_t bucketid1 =
                myTrie::hashRelative::hash(key, s->get_length(), 1) %
                Bucket_num;
            size_t bucketid2 =
                myTrie::hashRelative::hash(key, s->get_length(), 2) %
                Bucket_num;
            return current_bucketid == bucketid1 ? bucketid2 : bucketid1;
        }

        // Return a empty slot_id in bucketid
        int cuckoo_hash(size_t bucketid, htrie_map<CharT, T>* hm) {
            rehash_total_num++;
            uint64_t sta = get_time();

            // Set up the backup for recovery if the cuckoo hash fail
            slot* key_metas_backup = new slot[Bucket_num * cur_associativity]();
            memcpy(key_metas_backup, key_metas,
                   Bucket_num * cur_associativity * sizeof(slot));

            int ret_index = -1;

            /*
             * key_metas:   | x | x | x | x |   extra_slot: | y |
             *              | x | x | x | x |   slot wait to be exchange
             *              | x | x | x | x |
             */
            slot* extra_slot = new slot(0, 0, 0, 0);

            // cur_process_bucketid, cur_process_slotid indicate the
            // extra_slot's destination
            int cur_process_bucketid = bucketid;
            int cur_process_slotid = cur_associativity - 1;

            page_group_package pgp =
                hm->pm->get_page_group_package(normal_pgid, special_pgid);

            map<T, int> searchPoint_wait_2_be_update;
            for (int cuckoo_hash_time = 0; cuckoo_hash_time != Max_loop;
                 cuckoo_hash_time++) {
                /*
                 * The get_previous_slotid_in_same_bucket() will cause the -1 if
                 * there aren't any available slot can be cuckoo hash
                 * If the cur_process_slotid equals to -1, it means that all
                 * elements in current bucket's is anti-moved
                 */
                if (cur_process_slotid == -1) break;

                // Get the slot* we are replacing destination
                int cur_process_index =
                    get_index(cur_process_bucketid, cur_process_slotid);
                slot* cur_process_slot = get_slot(cur_process_index);

                /* Check that whether the cur_process_slot is anti-moved */
                // Get the another bucketid the cur_process_slot can be at
                int cur_kick_to_bucketid = get_another_bucketid(
                    pgp, cur_process_slot, cur_process_bucketid);
                // If the cur_process_bucketid == cur_kick_to_bucketid, we
                // process previous slotid
                if (cur_process_bucketid == cur_kick_to_bucketid ||
                    cur_process_index == ret_index) {
                    cur_process_slotid =
                        get_previous_slotid_in_same_bucket(cur_process_slotid);
                    continue;
                }

                /* Checking work is done, executing the cuckoo hashing */
                // Swap the extra slot content with the current process slot
                extra_slot->swap(cur_process_slot);

                // Add this value=index for the searchPoint index updateing
                searchPoint_wait_2_be_update[pgp.get_value(cur_process_slot)] =
                    get_column_store_index(cur_process_slot);

                // The first time swap the extra_slot indicate the
                // cur_process_slotid is ret_index
                if (ret_index == -1) {
                    ret_index = cur_process_index;
                }

                // cur_process_slot is a empty slot, cuckoo hash is done
                if (extra_slot->is_empty()) {
                    delete[] key_metas_backup;
                    delete extra_slot;

                    hm->apply_the_changed_searchPoint(
                        searchPoint_wait_2_be_update);

                    rehash_cost_time += get_time() - sta;

                    // return slot_id
                    return ret_index % cur_associativity;
                }

                // update the current bucketid and slotid which are the
                // replacing destination in next iteration
                cur_process_bucketid = cur_kick_to_bucketid;
                cur_process_slotid =
                    get_last_slotid_in_bucket(cur_kick_to_bucketid);
            }

            // Recover the key_metas
            memcpy(key_metas, key_metas_backup,
                   Bucket_num * cur_associativity * sizeof(slot));

            delete[] key_metas_backup;
            delete extra_slot;

            rehash_cost_time += get_time() - sta;

            // The cuckoo hash time exceeds Max_loop return -1 as slotid to
            // indicate cuckoo hash failed
            return -1;
        }

        /*------------------ 3. bursting function------------------*/
        inline unsigned int round_up_2_next_power_2(unsigned int x) {
            x--;
            x |= x >> 1;   // handle  2 bit numbers
            x |= x >> 2;   // handle  4 bit numbers
            x |= x >> 4;   // handle  8 bit numbers
            x |= x >> 8;   // handle 16 bit numbers
            x |= x >> 16;  // handle 32 bit numbers
            x++;
            return x;
        }

        inline unsigned int get_expected_associativity(unsigned int need_size) {
            return round_up_2_next_power_2(
                (double)need_size / (double)Bucket_num / Burst_ratio + 1);
        }

        inline bool need_burst() const {
            return elem_num >= Max_slot_num * Burst_ratio;
        }

        /*----------------searching in hash_node----------------*/
        void move_suffix_to_new_page(htrie_map<CharT, T>* hm,
                                  vector<page>& new_normal_page,
                                  vector<page>& new_special_page) {
            for (int i = 0; i != Bucket_num; i++) {
                for (int j = 0; j != cur_associativity; j++) {
                    slot* s = get_slot(i, j);

                    if (s->is_empty()) break;

                    if (s->is_special())
                        s->set_slot(hm->write_kv_to_page(
                            hm->get_content_pointer(s), s->get_length(),
                            hm->get_value(s), new_special_page));
                    else
                        s->set_slot(hm->write_kv_to_page(
                            hm->get_content_pointer(s), s->get_length(),
                            hm->get_value(s), new_normal_page));
                }
            }
        }

        /*----------------searching in hash_node----------------*/

        std::pair<bool, T> find_kv_in_pages(htrie_map<CharT, T>* hm, slot* s,
                                            const CharT* key,
                                            size_t keysize) {
            if (myTrie::hashRelative::keyEqual(
                    hm->get_content_pointer(get_page_group_id(s), s),
                                               s->get_length(), key, keysize)) {
                return std::pair<bool, T>(
                    true, *((T*)(hm->get_content_pointer(get_page_group_id(s), s) +
                                 s->get_length())));
            }
            return std::pair<bool, T>(false, T());
        }

        // return <found?, iterator>
        // iterator:    if slotid==-1, bucket is full
        //              if slotid!=-1, slotid is the insert position
        std::pair<bool, iterator> find_in_bucket(htrie_map<CharT, T>* hm,
                                                 size_t bucketid,
                                                 const CharT* key,
                                                 size_t keysize) {
            // print_key_metas(hm);
            // find the hitted slot in hashnode
            for (int i = 0; i != cur_associativity; i++) {
                slot* target_slot = get_slot(bucketid, i);
                if (target_slot->is_empty()) {
                    return std::pair<bool, iterator>(
                        false, iterator(false, T(), this, bucketid, i));
                }

                std::pair<bool, T> res =
                    find_kv_in_pages(hm, target_slot, key, keysize);
                if (res.first) {
                    return std::pair<bool, iterator>(
                        true, iterator(true, res.second, this, bucketid, i));
                }
            }
            return std::pair<bool, iterator>(
                false, iterator(false, T(), this, bucketid, -1));
        }

        iterator search_kv_in_hashnode(const CharT* key, size_t keysize,
                                       htrie_map<CharT, T>* hm) {
            // if found the existed target in bucket1 or bucket2, just
            // return the iterator for being modified or read
            size_t bucketId1 =
                myTrie::hashRelative::hash(key, keysize, 1) % Bucket_num;
            std::pair<bool, iterator> res1 =
                find_in_bucket(hm, bucketId1, key, keysize);

            if (res1.first) {
                return res1.second;
            }

            size_t bucketId2 =
                myTrie::hashRelative::hash(key, keysize, 2) % Bucket_num;
            std::pair<bool, iterator> res2 =
                find_in_bucket(hm, bucketId2, key, keysize);

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

        inline size_t get_page_group_id(size_t keysize) {
            return get_group_type(keysize) == group_type::SPECIAL_GROUP
                       ? special_pgid
                       : normal_pgid;
        }

        inline size_t get_page_group_id(slot* sl) {
            return get_group_type(sl->get_length()) == group_type::SPECIAL_GROUP
                       ? special_pgid
                       : normal_pgid;
        }

        void insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, size_t bucketid,
                                                 int slotid) {
            // if slotid==-1, it denotes that the bucket(bucketid) is full ,
            // so we rehash the key_metas
            if (slotid == -1) {
#ifdef REHASH_BEFORE_EXPAND
                if ((slotid = cuckoo_hash(bucketid, hm)) == -1) {
                    bool expand_success = expand_key_metas_space();
                    if (!expand_success) {
                        return std::pair<bool, T>(false, T());
                    } else {
                        // if expand success, we get new elem a empty slot
                        // in bucketid
                        for (int i = 0; i != cur_associativity; i++) {
                            slot* empty_slot =
                                key_metas + bucketid * cur_associativity + i;
                            if (empty_slot->is_empty()) {
                                slotid = i;
                                break;
                            }
                        }
                    }
                }
#else
                bool expand_success = expand_key_metas_space();
                if (!expand_success) {
                    if ((slotid = cuckoo_hash(bucketid, hm)) == -1) {

                        const string &prefix = this->anode::get_prefix();
                        trie_node* new_parent = hm->burst(burst_package(this, key_metas, Bucket_num, cur_associativity,
                                        hm->pm->get_page_group_package(
                                            normal_pgid, special_pgid)),
                            hm, this->anode::get_parent(),
                            prefix);

                        hm->access_kv_in_htrie_map(new_parent, key, keysize, v, false, prefix.data(), prefix.size());
                        delete this;
                        return;
                    }
                } else {
                    // if expand success, we get new elem a empty slot in
                    // bucketid
                    for (int i = 0; i != cur_associativity; i++) {
                        slot* empty_slot =
                            key_metas + bucketid * cur_associativity + i;
                        if (empty_slot->is_empty()) {
                            slotid = i;
                            break;
                        }
                    }
                    }
#endif
            }

            // now the slotid cannot be -1 and slotid is lower than
            // Associativity
            assert(slotid != -1 && slotid >= 0 && slotid < cur_associativity);

            slot* target_slot = get_slot(bucketid, slotid);

            // Ask page_manager that whether there is a place for new element
            /*
             * if page_manager pm don't have enough memory, pm will take charge
             * of the page group allocation, and update all hash_nodes'
             * normal_pgid, special_pgid. And so the hm->write_kv() can write
             * the content to a updated page group*/
            /*
             * if page_manager pm have enough memory, the hm->write_kv() will
             * write the content to original page group
             */
            hm->pm->try_insert(get_group_type(keysize),
                               get_page_group_id(keysize), keysize, hm);

            // call htrie-map function: write_kv_to_page ()
            // return a slot with position that element been written
            target_slot->set_slot(
                hm->write_kv(get_page_group_id(keysize), key, keysize, v));

            // set v2k
            hm->set_v2k(v, this, get_column_store_index(target_slot));

            elem_num++;
            
            return;
        }
    };

    /*
    * burst()
    * when the hash_node's element reach the need_busrt() threshold, we
    * burst this element into a tree with several small size of
    * hash_node or trie_node linking with those hash_node
    *
    * hash_node(size:100)
    *
    *       burst()↓
    *
    * trie_node(size:0)
    * childs: |hash_node(size:25)|trie_node(size:1)|hash_node(size:25)|)
    *                              |
    *                      hash_node(size:49)
    */
    struct burst_package {
       private:
        hash_node* bursting_node_;
        page_group_package pgp_;
        vector<slot> elems_;

       public:
        burst_package(hash_node* bursting_node, slot* elems, size_t bucket_num,
                      size_t associativity, page_group_package pgp)
            : bursting_node_(bursting_node), pgp_(pgp) {
          for (int i = 0; i != bucket_num; i++)
            for (int j = 0; j != associativity; j++) {
              slot* s = elems + i * associativity + j;
              if (s->is_empty()) break;
              elems_.push_back(*s);
            }
        }

        void update_burst_package(page_manager* new_pm, group_type resize_type) {
            size_t n_group_id = new_pm->require_group_id(group_type::NORMAL_GROUP);
            size_t s_group_id =new_pm->require_group_id(group_type::SPECIAL_GROUP);

            page_group_package new_pgp = new_pm->get_page_group_package(n_group_id, s_group_id);
            for (auto& s : elems_) {
                // ignore the slot that not belong to current resize group
                // type
                if (get_group_type(&s) != resize_type) continue;

                // get the content from old page_manager and write it to the
                // new page_manager
                s.set_slot(new_pm->write_kv((resize_type == group_type::NORMAL_GROUP
                                        ? n_group_id
                                        : s_group_id),
                                    pgp_.get_content_pointer(&s),
                                    s.get_length(), pgp_.get_value(&s)));
            }
            pgp_.set_page_group(resize_type, new_pgp.get_page_group(resize_type));
        }

        slot operator[](int index) { return elems_[index]; }
        const slot top() const { return elems_.back(); }
        void pop() { elems_.pop_back(); }

        size_t size() {return elems_.size(); }
        page_group_package get_pgp() const { return pgp_; }
        const hash_node* get_bursting_node() const {return bursting_node_; }

        void print_bp() {
            cout << "-----------------\n";
            for(int i=0; i!= elems_.size(); i++){
                elems_[i].print_slot(pgp_);
            }
            cout << endl;
        }
        
        inline unsigned int cal_common_prefix_len(const char* s1,
                                           unsigned int cur_longest_prefix_len,
                                           const char* s2,
                                           unsigned int new_key_size) {
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

        string get_common_prefix() {
            // calculate the capacity of hashnode we need
            const char* ret_key_pointer = nullptr;
            unsigned int common_prefix_len = INT_MAX;
            for (int i = 0; i != elems_.size() && common_prefix_len != 0; i++){
                    char* key = pgp_.get_content_pointer(&(elems_[i]));
                    if (ret_key_pointer == nullptr) ret_key_pointer = key;

                    // update the common_prefix_len
                    unsigned int cur_com_prefix_len = cal_common_prefix_len(
                        ret_key_pointer, common_prefix_len, key,
                        elems_[i].get_length());
                    if (common_prefix_len > cur_com_prefix_len) {
                        common_prefix_len = cur_com_prefix_len;
                }
            }
            return string(ret_key_pointer, common_prefix_len);
        }
    };

    // bursting function
    trie_node* burst(burst_package bp, htrie_map* hm, trie_node* orig_parent, const string &orig_prefix) {
        burst_total_counter++;

        pm->register_burst_package(&bp);

        // The return header
        trie_node* ret_trie_root = new trie_node(orig_parent, orig_prefix.data(), orig_prefix.size());

        if (orig_parent == nullptr)
            hm->set_t_root(ret_trie_root);
        else
            orig_parent->add_child(orig_prefix.back(), ret_trie_root);

        trie_node* parent = ret_trie_root;

        // Get the common_prefix to eliminate the redundant burst
        string common_prefix = bp.get_common_prefix();
        const char* common_prefix_key = common_prefix.data();
        unsigned int common_prefix_keysize = common_prefix.size();

        // New prefix = prior prefix + common chain prefix
        string prefix = orig_prefix + common_prefix;

        // Create the common prefix trie chain with several single trie_node
        // The number of node is common_prefix_keysize
        for (int i = 0; i != common_prefix_keysize; i++) {
            trie_node* cur_trie_node =
                new trie_node(parent, prefix.data(),
                                prefix.size() - common_prefix_keysize + i + 1);
            parent->add_child(common_prefix_key[i], cur_trie_node);
            parent = cur_trie_node;
        }

        // Insert the elements with same first char after common_prefix_len
        while (bp.size() != 0) {
                slot s = bp.top();

                char* new_key = bp.get_pgp().get_content_pointer(&s) + common_prefix_keysize;
                size_t length_left = s.get_length() - common_prefix_keysize;
                T v = bp.get_pgp().get_value(&s);

                hm->access_kv_in_htrie_map(parent, new_key, length_left, v, false, prefix.data(), prefix.size());

                bp.pop();
        }

        // remove current hash_node in the notify list in page_manager
        pm->remove_burst_package(&bp);

        return ret_trie_root;
    }

    class page_manager {
       public:
        class page_group {
            class page {
               public:
                unsigned int cur_pos;
                char* content;

                page() : cur_pos(0), content(nullptr) {}

                void init_page(size_t size_per_page) {
                    content = (char*)malloc(size_per_page);
                }

                void append_impl(const CharT* key, size_t keysize, T& value,
                                 unsigned int alignment = 1) {
                    // append the string
                    std::memcpy(content + cur_pos, key,
                                keysize * sizeof(CharT));

                    // append the value
                    std::memcpy(content + cur_pos + keysize * sizeof(CharT),
                                &value, sizeof(T));

                    cur_pos += calc_align(keysize * sizeof(CharT) + sizeof(T),
                                          alignment);
                }

                ~page() {
                    if (content != nullptr) {
                        free(content);
                        // static int counting = 0;
                        // cout << "page dec " << counting++ << endl;
                    }
                }
            };

            page* pages;
            int cur_page_id;
            bool is_special;

           public:
            page_group() : pages(nullptr), cur_page_id(-1), is_special(false) {}

            void init_pg(int page_number, bool spe) {
                is_special = spe;
                cur_page_id = 0;
                pages = new page[spe ? MAX_SPECIAL_PAGE : MAX_NORMAL_PAGE]();
                pages[0].init_page(spe ? DEFAULT_SPECIAL_Max_bytes_per_kv
                                       : DEFAULT_Max_bytes_per_kv);
            }

            // get function
            inline char* get_content_pointer_in_page(slot* s) {
                return pages[s->get_page_id()].content + s->get_pos();
            }

            inline T get_value_in_page(slot* s) {
                T v;
                std::memcpy(&v, get_content_pointer_in_page(s) + s->get_length(),
                            sizeof(T));
                return v;
            }

            // set function
            slot write_kv_to_page(const CharT* key, size_t keysize, T v) {
                // allocate space
                size_t need_size = keysize * sizeof(CharT) + sizeof(T);

                if (pages[cur_page_id].cur_pos + need_size >
                    (is_special ? DEFAULT_SPECIAL_Max_bytes_per_kv
                                : Max_bytes_per_kv)) {
                    cur_page_id++;
                    pages[cur_page_id].init_page(
                        is_special ? DEFAULT_SPECIAL_Max_bytes_per_kv
                                   : DEFAULT_Max_bytes_per_kv);
                }

                // get page being written
                page& target_page = pages[cur_page_id];

                // record position before updating and status modify
                slot ret_slot =
                    slot(is_special, keysize, target_page.cur_pos, cur_page_id);

                // write content
                target_page.append_impl(key, keysize, v,
                                        is_special ? DEFAULT_SPECIAL_ALIGNMENT
                                                   : DEFAULT_NORMAL_ALIGNMENT);

                return ret_slot;
            }

            page_group(const page_group& orig) {
                cout << "calling copy func\n";
                pages = orig.pages;
                cur_page_id = orig.cur_page_id;
                is_special = orig.is_special;
            }

            inline size_t get_cur_page_id(){
                return cur_page_id;
            }

            inline size_t get_max_page_id(){
                return is_special ? MAX_SPECIAL_PAGE : MAX_NORMAL_PAGE;
            }

            inline size_t get_max_per_page_size() {
                return is_special ? DEFAULT_SPECIAL_Max_bytes_per_kv : DEFAULT_Max_bytes_per_kv;
            }

            bool try_insert(size_t try_insert_keysize) {
                if (cur_page_id + 1 < get_max_page_id()) return true;
                if ((pages[cur_page_id].cur_pos +
                     try_insert_keysize * sizeof(CharT) + sizeof(T)) <=
                    get_max_per_page_size())
                    return true;
                return false;
            }

            ~page_group() {
                // cout << (void*)this << "page_group deconstructor" << endl;
                int pages_size =
                    is_special ? MAX_SPECIAL_PAGE : MAX_NORMAL_PAGE;
                delete []pages;
            }
        };

        page_group* normal_pg;
        page_group* special_pg;

        size_t n_size;
        size_t s_size;

        void init_a_new_page_group(group_type type, size_t page_group_index) {
            if (type == group_type::SPECIAL_GROUP) {
                s_size++;
                special_pg[page_group_index].init_pg(MAX_SPECIAL_PAGE, true);
                // cout << "special page group increase!" << endl;
                return;
            } else if (type == group_type::NORMAL_GROUP) {
                n_size++;
                normal_pg[page_group_index].init_pg(MAX_NORMAL_PAGE, false);
                // cout << "normal page group increase!" << endl;
                return;
            } else {
                cout << "undefined type!" << endl;
                assert(false);
                exit(0);
                return;
            }
        }

       public:
        page_manager(size_t normal_page_group_number = 1,
                     size_t special_page_group_number = 1)
            : normal_pg(new page_group[128]),
              special_pg(new page_group[128]),
              n_size(0),
              s_size(0) {
            for (int i = 0; i != normal_page_group_number; i++)
                init_a_new_page_group(group_type::NORMAL_GROUP, i);

            for (int i = 0; i != special_page_group_number; i++)
                init_a_new_page_group(group_type::SPECIAL_GROUP, i);
        }

        ~page_manager(){
            delete []normal_pg;
            delete []special_pg;
            // cout << "page manager deconstruction" << endl;
        }

        inline slot write_kv(size_t page_group_id, const CharT* key,
                             size_t keysize, T v) {
            return get_group_type(keysize) == group_type::SPECIAL_GROUP
                       ? special_pg[page_group_id].write_kv_to_page(key,
                                                                    keysize, v)
                       : normal_pg[page_group_id].write_kv_to_page(key, keysize,
                                                                   v);
        }

        inline size_t require_group_id(group_type gt) {
            size_t least_page_page_group_id = 0;
            size_t least_page = SIZE_MAX;
            page_group* pgs =
                (gt == group_type::NORMAL_GROUP ? normal_pg : special_pg);
            for (int i = 0; i != n_size; i++) {
                size_t cur_least_page = pgs[i].get_cur_page_id();
                if (least_page > cur_least_page) {
                    least_page = cur_least_page;
                    least_page_page_group_id = i;
                }
            }
            return least_page_page_group_id;
        }

        inline char* get_content_pointer_in_pm(size_t page_group_id,
                                            slot* s) {
            return get_group_type(s) == group_type::SPECIAL_GROUP
                       ? special_pg[page_group_id].get_content_pointer_in_page(s)
                       : normal_pg[page_group_id].get_content_pointer_in_page(s);
        }

        inline T get_value_in_pm(size_t page_group_id, slot* s) {
            return get_group_type(s) == group_type::SPECIAL_GROUP
                       ? special_pg[page_group_id].get_value_in_page(s)
                       : normal_pg[page_group_id].get_value_in_page(s);
        }

        inline page_group_package get_page_group_package(size_t n_pg,
                                                         size_t s_pg) {
            return page_group_package(normal_pg + n_pg, special_pg + s_pg);
        }

        void set_n_size(size_t new_n_size) { n_size = new_n_size; }

        void set_s_size(size_t new_s_size) { s_size = new_s_size; }

        void set_normal_pg(page_group* temp_normal_pg) {
            normal_pg = temp_normal_pg;
        }

        void set_special_pg(page_group* temp_special_pg) {
            special_pg = temp_special_pg;
        }

        void swap(page_manager* pm, group_type gt) {
            // swap the normal part
            // swap the normal page group
            if (gt == group_type::NORMAL_GROUP) {
                page_group* temp_normal_pg = normal_pg;
                normal_pg = pm->normal_pg;
                pm->set_normal_pg(temp_normal_pg);

                // swap the n_size
                int temp_n_size = n_size;
                set_n_size(pm->n_size);
                pm->set_n_size(temp_n_size);
                return;
            }

            // swap the special part
            // swap the special page group
            page_group* temp_special_pg = special_pg;
            special_pg = pm->special_pg;
            pm->set_special_pg(temp_special_pg);

            // swap the s_size
            int temp_s_size = s_size;
            set_s_size(pm->s_size);
            pm->set_s_size(temp_s_size);
            return;
        }

        void swap(page_manager &pm){
            swap(pm, group_type::NORMAL_GROUP);
            swap(pm, group_type::SPECIAL_GROUP);
        }

        /* 
         * Oberserver design pattern
         * page_manager is a Subject that have the function of register, remove,
         * notify(traverse_for_pgm_resize) those burst_package in burst()
         */
        vector<burst_package*> notify_list;
        void register_burst_package(burst_package *add_bp_ptr) {
            notify_list.push_back(add_bp_ptr);
        }

        void remove_burst_package(const burst_package *const rm_bp_ptr){
            for(auto it = notify_list.begin(); it!=notify_list.end(); it++) {
                if( (*it)->get_bursting_node() == rm_bp_ptr->get_bursting_node()){ 
                    notify_list.erase(it);
                    break;
                }
            }
        }

        void notify_burst_package(page_manager *new_pm, group_type resize_type) {
            for (auto bp_ptr : notify_list) {
                bp_ptr->update_burst_package(new_pm, resize_type);
            }
        }

        void clean_useless_in_pm(group_type resize_type,
                                 htrie_map<CharT, T>* hm,
                                 size_t expand_ratio = 1) {
            uint64_t sta = get_time();

            page_manager* pm;
            // cout << "resizing!!!\n";
            if (resize_type == group_type::SPECIAL_GROUP) {
                pm = new page_manager(0, s_size << expand_ratio);
            } else {
                pm = new page_manager(n_size << expand_ratio, 0);
            }

            // try insert, if failed, we reallocate the page groups,
            // update the pgid in hashnodes and return
            anode* root = hm->t_root;
            root->traverse_for_pgm_resize(this, pm, resize_type);

            // notify the bursting hash node that your keymetas have been change
            // because of the resize
            notify_burst_package(pm, resize_type);

            // old page_manager <=swap=> new page_manager
            swap(pm, resize_type);

            uint64_t end = get_time();
            // cout << "resizeing cost: " << (end - sta) / 1000 / (double)1000
            //      << " s" << endl;

            delete pm;
            return;
        }

        void try_insert(group_type resize_type, size_t page_group_id,
                        size_t try_insert_keysize, htrie_map<CharT, T>* hm) {

            // get target page_group
            page_group& target = (resize_type == group_type::SPECIAL_GROUP
                                      ? special_pg[page_group_id]
                                      : normal_pg[page_group_id]);

            // try insert, if success just return
            if (target.try_insert(try_insert_keysize)) return;

            // traverse the trie and clean the useless content
            clean_useless_in_pm(resize_type, hm);
            return;
        }
    };

    inline size_t get_normal_group_id() {
        return pm->require_group_id(group_type::NORMAL_GROUP);
    }
    inline size_t get_special_group_id() {
        return pm->require_group_id(group_type::SPECIAL_GROUP);
    }

    inline char* get_content_pointer(size_t page_group_id, slot* s) {
        return pm->get_content_pointer_in_pm(page_group_id, s);
    }

    inline T get_value(size_t page_group_id, slot* s) {
        return pm->get_value_in_pm(page_group_id, s);
    }

    inline slot write_kv(size_t page_group_id, const CharT* key, size_t keysize,
                         T v) {
        return pm->write_kv(page_group_id, key, keysize, v);
    }

    class SearchPoint {
       public:
        anode* node;
        int index;

        SearchPoint() : node(nullptr), index(-1) {}
        SearchPoint(anode* n, int i) : node(n), index(i) {}

        void set_index(int i) { index = i; }

        std::string get_string(page_manager* pm) {
            if (node == nullptr) return string();
            // if the node is trie_node, just return the prefix on node
            if (node->is_trie_node()) {
                return node->get_prefix();
            }
            string res = node->get_prefix();
            if (index != -1) {
                hash_node* hnode = (hash_node*)node;
                slot* sl = hnode->get_column_store_slot(index);
                res = res + string(pm->get_content_pointer_in_pm(
                                       hnode->get_page_group_id(sl), sl),
                                   sl->get_length());
            }
            return res;
        }
    };

   public:

    void set_searchPoint_index(T v, int index) { v2k[v].set_index(index); }

    void set_v2k(T v, anode* node, int index) {
        v2k[v] = SearchPoint(node, index);
    }

    // function for batch updating the searchPoints to v2k
    void apply_the_changed_searchPoint(map<T, int>& searchPoints) {
        for (auto it = searchPoints.begin(); it != searchPoints.end(); it++)
            set_searchPoint_index(it->first, it->second);
    }

    // std::map<T, SearchPoint> v2k;
    boost::unordered_map<T, SearchPoint> v2k;

    anode* t_root;
    
    page_manager *pm;

    htrie_map(size_t customized_associativity = DEFAULT_Associativity,
              size_t customized_bucket_count = DEFAULT_Bucket_num,
              size_t customized_byte_per_kv = DEFAULT_Max_bytes_per_kv,
              double customized_burst_ratio = DEFAULT_Burst_ratio)
        : t_root(nullptr), pm(new page_manager()) {
        std::cout << "SET UP GROWING-CUCKOOHASH-TRIE MAP\n";
        cout << "GROW_ASSOCIATIVITY\n";
        cout << "PM\n";

#ifdef REHASH_BEFORE_EXPAND
        std::cout << "REHASH_BEFORE_EXPAND\n";

#else
        std::cout << "EXPAND_BEFORE_REHASH\n";

#endif

#ifdef MAP_REP
        cout << "MAP REPRESENTATION\n";
#endif
#ifdef ARRAY_REP
        cout << "ARRAY REPRESENTATION\n";
#endif
#ifdef LIST_REP
        cout << "LIST REPRESENTATION\n";
#endif

        Associativity = customized_associativity;
        Bucket_num = customized_bucket_count;
        Max_bytes_per_kv = customized_byte_per_kv;
        Burst_ratio = customized_burst_ratio;

        Max_slot_num = Associativity * Bucket_num;
        Max_loop = Max_slot_num * 0.5;

        t_root = new hash_node(nullptr, string(), this, Associativity);
    }

    // TODO deconstructor

    void set_t_root(anode* node) { t_root = node; }

    class slot {
       public:
        uint32_t encode;

        slot() : encode(0) {}

        // encode slot as | special | length | position | page_id |
        uint64_t encode_slot(bool spe, uint64_t len, uint64_t po, uint64_t pd) {
            encode = pd;
            if (spe) {
                assert(len < 1 << NBITS_LEN_S &&
                       (po / DEFAULT_SPECIAL_ALIGNMENT) < 1 << NBITS_POS_S &&
                       pd < 1 << NBITS_PID_S);
                encode += ((po / DEFAULT_SPECIAL_ALIGNMENT) << NBITS_PID_S);
                encode += len << (NBITS_PID_S + NBITS_POS_S);
                encode += ((uint64_t)1)
                          << (NBITS_PID_S + NBITS_LEN_S + NBITS_POS_S);
            } else {
                assert(len < 1 << NBITS_LEN &&
                       (po / DEFAULT_NORMAL_ALIGNMENT) < 1 << NBITS_POS &&
                       pd < 1 << NBITS_PID);
                // encode into a 64bit data
                encode += (po / DEFAULT_NORMAL_ALIGNMENT) << NBITS_PID;
                encode += len << (NBITS_PID + NBITS_POS);
                encode += ((uint64_t)0)  << (NBITS_PID + NBITS_LEN + NBITS_POS);
            }

            return encode;
        }

        bool is_empty() { return get_length() == 0; }

        bool is_special() { return get_special(); }

        slot(bool is_special, KeySizeT l, size_t p, size_t pi) {
            encode = encode_slot(is_special, l, p, pi);
        }

        void set_slot(bool is_special, KeySizeT l, size_t p, size_t pi) {
            encode = encode_slot(is_special, l, p, pi);
        }

        void set_slot(slot s) { encode = s.encode; }

        void set_slot(slot* s) { encode = s->encode; }

        bool get_special() {
            return (encode >> (NBITS_PID + NBITS_LEN + NBITS_POS)) == 1;
        }

        size_t get_length() {
            return encode >> (NBITS_LEN + NBITS_POS + NBITS_PID)
                       ? ((encode >> (NBITS_PID_S + NBITS_POS_S)) %
                          (1 << NBITS_LEN_S))
                       : ((encode >> (NBITS_PID + NBITS_POS)) %
                          (1 << NBITS_LEN));
        }

        size_t get_pos() {
            return encode >> (NBITS_PID + NBITS_LEN + NBITS_POS)
                       ? ((encode >> NBITS_PID_S) % (1 << NBITS_POS_S) *
                          DEFAULT_SPECIAL_ALIGNMENT)
                       : ((encode >> NBITS_PID) % (1 << NBITS_POS) *
                          DEFAULT_NORMAL_ALIGNMENT);
        }

        size_t get_page_id() {
            return encode >> (NBITS_PID + NBITS_LEN + NBITS_POS)
                       ? encode % (1 << NBITS_PID_S)
                       : encode % (1 << NBITS_PID);
        }

        void swap(slot *sl) {
            uint32_t temp_encode = sl->encode;
            sl->encode = encode;
            encode = temp_encode;
        }

        void print_slot(page_group_package &pgp) {
            cout << get_special() << "," << get_length() << "," << get_pos()
                << "," << get_page_id() << ","
                << string(pgp.get_content_pointer(this), get_length()) << ","
                << pgp.get_value(this) << endl;
        }
    };

    inline static unsigned int calc_align(unsigned int n, unsigned align) {
        return ((n + align - 1) & (~(align - 1)));
    }

    class iterator {
       public:
        bool found;
        const T v;
        anode* target_node;
        size_t bucketid;
        int slotid;

        iterator(bool f, T vv, anode* hnode, size_t bid, int sid)
            : found(f), v(vv), target_node(hnode), bucketid(bid), slotid(sid) {}

        void insert_hashnode(const CharT* key, size_t key_size,
                                           htrie_map<CharT, T>* hm, T v) {
            ((hash_node*)target_node)
                ->insert_kv_in_hashnode(key, key_size, hm, v, bucketid, slotid);
            return;
        }
    };

    // TODO: divided into find and insert mode
    std::pair<bool, T> access_kv_in_htrie_map(anode* start_node,
                                              const CharT* key, size_t key_size,
                                              T v, bool findMode,
                                              const CharT* prefix_key = nullptr,
                                              size_t prefix_key_size = 0) {
      // update longest_string_size
      longest_string_size =
          longest_string_size > key_size ? longest_string_size : key_size;

      anode* current_node = start_node;

      // TODO: pos updating need refine?
      // The pos update is moved to find_trie_node_child(fast-path or
      // non-fast-path way) while the pos increment
      for (size_t ref_pos = 0; ref_pos < key_size;) {
        switch (current_node->get_node_type()) {
          case node_type::TRIE_NODE: {
              trie_node* orig_tnode = (trie_node*)current_node;
            // return the hitted trie_node* or create a new
            // trie_node with a hash_node son
            current_node = ((trie_node*)current_node)
                               ->find_trie_node_child(findMode, key, ref_pos,
                                                      key_size, this);

            if(current_node == nullptr){
              if (findMode) {
                return std::pair<bool, T>(false, T());
              } else {
                string new_prefix = string(prefix_key, prefix_key_size) + string(key, ref_pos);
                // Create a corresponding hash_node and add it to current
                // trie_node's child representation
                current_node =
                    new hash_node(orig_tnode, new_prefix, this);
                orig_tnode->add_child(key[ref_pos - 1], current_node);
              }
            } 

          } break;
          case node_type::HASH_NODE: {
            iterator it = ((hash_node*)current_node)
                              ->search_kv_in_hashnode(key + ref_pos,
                                                      key_size - ref_pos, this);
            if (findMode) {
              return std::pair<bool, T>(it.found, it.v);
            } else {
              it.insert_hashnode(key + ref_pos, key_size - ref_pos, this, v);
              return std::pair<bool, T>(it.found, it.v);
            }
          } break;
          default:
            cout << "wrong type!";
            exit(0);
        }
      }

      // find a key in node's only value
      iterator it = current_node->search_kv_in_node();

      if (findMode) {
        return std::pair<bool, T>(it.found, it.v);
      } else {
        current_node->insert_value_in_node(string(prefix_key, prefix_key_size) + string(key, key_size), v, this);
        return std::pair<bool, T>(true, v);
      }
    }

    /*---------------external accessing interface-------------------*/

    // search operation
    T searchByKey(std::string key) {
        return access_kv_in_htrie_map(t_root, key.data(), key.size(), T(), true).second;
    }

    std::string searchByValue(T v) { return v2k[v].get_string(pm); }

    // find operation
    std::pair<bool, T> findByKey(std::string key) {
        return access_kv_in_htrie_map(t_root, key.data(), key.size(), T(), true);
    }

    std::pair<bool, std::string> findByValue(T v) {
        if (v2k.find(v) == v2k.end()) {
            return std::pair<bool, T>(false, string());
        } else {
            return std::pair<bool, T>(true, v2k[v].get_string(this));
        }
    }

    // insert operation
    std::pair<bool, T> insertKV(std::string key, T v) {
        return access_kv_in_htrie_map(t_root, key.data(), key.size(), v, false);
    }

    /*---------------external cleaning interface-------------------*/
    void clean_useless() {
        // zero at last means that we don't need to expand the page_manager
        pm->clean_useless_in_pm(group_type ::NORMAL_GROUP, this, 0);
        pm->clean_useless_in_pm(group_type ::SPECIAL_GROUP, this, 0);
    }

    void traverse_level() {
        cout << endl;
        queue<anode*> q;
        q.push(t_root);
        unsigned int fpm_memory = 0;
        while (!q.empty()) {
            anode* cur_node = q.front();
            q.pop();
            if (cur_node->is_hash_node()) continue;

            unsigned int cur_fpm_mem =
                ((trie_node *)cur_node)->fpm_ == nullptr
                    ? 0
                    : ((trie_node *)cur_node)->fpm_->get_fpm_memory();
            fpm_memory += cur_fpm_mem;
            
            vector<anode*> childs;
            ((trie_node*)cur_node)->get_childs_vector(childs);
            for (auto c : childs) q.push(c);
        }
        cout << "fast path memory: " << fpm_memory << endl;
    }
};  // namespace myTrie

}  // namespace myTrie