#include "../util/my_timer.hpp"

#include <iostream>
#include <vector>
#include <queue>
#include <string>

#include <assert.h>

#include <boost/unordered_map.hpp>

/*---- Default configuration ---*/
static const unsigned int DEFAULT_ASSOCIATIVITY = 8;
static const unsigned int DEFAULT_BUCKET_NUM = 59;
static const unsigned int DEFAULT_NORMAL_PAGE_SIZE = 4096;
static const unsigned int DEFAULT_SPECIAL_PAGE_SIZE = (4096 * 4);
static const double DEFAULT_CUCKOO_HASH_RATIO = 0.5;
/*---- slot use bits configuration ---*/
enum { NBITS_SPECIAL = 1, NBITS_SPECIAL_S = 1 };  // is special
enum { NBITS_LEN = 7, NBITS_LEN_S = 13 };         // length
enum { NBITS_PID = 12, NBITS_PID_S = 8 };         // page id
enum { NBITS_POS = 12, NBITS_POS_S = 9 };         // position in page
// normal/special bound
static const unsigned int MAX_NORMAL_LEN = (1 << NBITS_LEN);
/*---- page manager configuration ---*/
static const unsigned int DEFAULT_NORMAL_PAGE_NUMBER = (1 << NBITS_PID);
static const unsigned int DEFAULT_SPECIAL_PAGE_NUMBER = (1 << NBITS_PID_S);
static const unsigned int DEFAULT_SPECIAL_ALIGNMENT = 32;
static const unsigned int DEFAULT_NORMAL_ALIGNMENT = 1;
/*---- fast path configuration ---*/
static const unsigned int FAST_PATH_NODE_NUM = 20;

/*---- For information ---*/
// todo: wait to be deleted, just for recording the time that expand() cost
// burst
uint32_t burst_total_counter = 0;
uint64_t burst_total_time = 0;
// expand
uint64_t expand_cost_time = 0;
// cuckoo hash
uint64_t cuckoohash_cost_time = 0;
uint64_t cuckoohash_total_num = 0;
// fast path establish
uint64_t establish_fastpath_cost_time = 0;

// TODO: myTrie, htrie_map rename: 
// myTrie: zyz_trie, htrie_map: bi_trie
namespace myTrie {
using namespace std;
// charT = char, T = value type
template <class CharT, class T, size_t BUCKET_NUM = DEFAULT_BUCKET_NUM, size_t ASSOCIATIVITY = DEFAULT_ASSOCIATIVITY>
class htrie_map {
   private:
    static bool key_equal(const CharT* key_lhs, std::size_t key_size_lhs,
                const CharT* key_rhs, std::size_t key_size_rhs) {
        if (key_size_lhs == 0 && key_size_rhs == 0) {
            return true;
        }
        if (key_size_lhs != key_size_rhs) {
            return false;
        } else {
            return std::memcmp(key_lhs, key_rhs, key_size_lhs * sizeof(CharT)) ==
                0;
        }
    }

    // For fasthash64:
    static inline uint64_t mix(uint64_t h) {
        h ^= h >> 23;
        h *= 0x2127599bf4325c37ULL;
        h ^= h >> 47;
        return h;
    }

    // A default hash function:
    static uint64_t fasthash64(const char* buf, size_t len, uint64_t seed) {
        uint64_t const m = 0x880355f21e6d1965ULL;
        uint64_t const* pos = (uint64_t const*)buf;
        uint64_t const* end = pos + (len / 8);
        const unsigned char* pos2;
        uint64_t h = seed ^ (len * m);
        uint64_t v;

        while (pos != end) {
            v = *pos++;
            h ^= mix(v);
            h *= m;
        }

        pos2 = (const unsigned char*)pos;
        v = 0;

        switch (len & 7) {
            case 7:
                v ^= (uint64_t)pos2[6] << 48;
            case 6:
                v ^= (uint64_t)pos2[5] << 40;
            case 5:
                v ^= (uint64_t)pos2[4] << 32;
            case 4:
                v ^= (uint64_t)pos2[3] << 24;
            case 3:
                v ^= (uint64_t)pos2[2] << 16;
            case 2:
                v ^= (uint64_t)pos2[1] << 8;
            case 1:
                v ^= (uint64_t)pos2[0];
                h ^= mix(v);
                h *= m;
        }

        return mix(h);
    }

    static size_t hash(const CharT* key, std::size_t key_size, size_t hashType = 1) {
        uint64_t hash;
        if (hashType == 1) {
            hash = fasthash64(key, key_size, 0xdeadbeefdeadbeefULL);
            return hash;
        } else {
            hash = fasthash64(key, key_size, 0xabcdefabcdef1234ULL);
            return hash;
        }
    }

    class trie_node;
    class hash_node;

    class found_result;

    class slot;

    /* page management */
    class page_manager;
    class page_manager_agent;
    class page_group;
    class page;

    /* group type divided by the keysize */
    enum class group_type : unsigned char { NORMAL_GROUP, SPECIAL_GROUP };

    /* node type definition */
    enum class node_type : unsigned char { HASH_NODE, TRIE_NODE };

    /* helper function */
    static inline group_type get_group_type(size_t keysize) {
        return keysize < MAX_NORMAL_LEN ? group_type::NORMAL_GROUP
                                        : group_type::SPECIAL_GROUP;
    }

    static inline group_type get_group_type(slot* s) {
        return s->get_length() < MAX_NORMAL_LEN ? group_type::NORMAL_GROUP
                                                : group_type::SPECIAL_GROUP;
    }

   private:
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

        slot(bool is_special, size_t l, size_t p, size_t pi) {
            encode = encode_slot(is_special, l, p, pi);
        }

        void set_slot(bool is_special, size_t l, size_t p, size_t pi) {
            encode = encode_slot(is_special, l, p, pi);
        }

        void set_slot(slot s) { encode = s.encode; }

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

        void print_slot(page_manager_agent &pm_agent) {
            cout << get_special() << "," << get_length() << "," << get_pos()
                << "," << get_page_id() << ","
                << string(pm_agent.get_content_pointer(this), get_length()) << ","
                << pm_agent.get_value(this) << endl;
        }
    };

    /* helper class for hash_node get its page_groups */
    class page_manager_agent {
        typename page_manager::page_group* n_group;
        typename page_manager::page_group* s_group;

       public:
        page_manager_agent(typename page_manager::page_group* ng,
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

        // Try to insert element to its right group and return availability
        inline bool try_insert(size_t keysize) {
            return get_group_type(keysize) == group_type::SPECIAL_GROUP
                       ? s_group->try_insert(keysize)
                       : n_group->try_insert(keysize);
        }

        // Insert element to its right group and return the slot(position)
        inline slot insert_element(const CharT* key, size_t keysize, T v) {
            return get_group_type(keysize) == group_type::SPECIAL_GROUP
                       ? s_group->write_kv_to_page(key, keysize, v)
                       : n_group->write_kv_to_page(key, keysize, v);
        }
    };

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

        found_result insert_value_in_node(const string &prefix, T v,
                                      htrie_map<CharT, T>* hm) {
            value_ = v;
            have_value_ = true;
            hm->set_v2k(v, this, -1);
            prefix_ = prefix;
            return found_result(have_value_, value_, -1, -1);
        }

        found_result search_kv_in_node() {
            return found_result(have_value_, value_, -1, -1);
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


    class trie_node : public anode {
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
                fast_path new_fast_path(hash(key, key_size, 1), string(key, key_size), node_ptr);

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
                unsigned int target_hash_val = hash(key, key_size);

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
                    if (key_equal(
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

        /*-----------------list child_representation-----------------*/
        class child_representation {
        private:
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

        ~trie_node() { if (fpm_ != nullptr) delete fpm_; }

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
        void print_slot(int i, int j, page_manager_agent &pm_agent) {
            slot* s = key_metas + i * cur_associativity + j;
            cout << i * cur_associativity + j << ":" << s->get_special() << ","
                 << s->get_length() << "," << s->get_pos() << ","
                 << s->get_page_id() << ",";
            string str = string(pm_agent.get_content_pointer(s),
                                s->get_length());
            cout << str;
            T v = pm_agent.get_value(s);
            cout << "=" << v << "\n";
        }

        static void print_slot(slot* s) {
            cout << s->get_special() << "," << s->get_length() << ","
                 << s->get_pos() << "," << s->get_page_id() << "," << endl;
        }

        void print_key_metas(htrie_map<CharT, T>* hm) {
            page_manager_agent pm_agent = hm->pm->get_page_manager_agent(normal_pgid, special_pgid);
            for (int i = 0; i != BUCKET_NUM; i++) {
                for (int j = 0; j != cur_associativity; j++) {
                    print_slot(i, j, pm_agent);
                }
                cout << "---\n";
            }
            cout << endl;
        }

       public:
        explicit hash_node(trie_node* p,const string &prefix, page_manager* pm,
                           size_t need_associativity = 1)
            : anode(node_type::HASH_NODE, p, prefix.data(), prefix.size()),
              cur_associativity(need_associativity > ASSOCIATIVITY
                                    ? ASSOCIATIVITY
                                    : need_associativity),
              elem_num(0),
              normal_pgid(pm->require_group_id(group_type::NORMAL_GROUP)),
              special_pgid(pm->require_group_id(group_type::SPECIAL_GROUP)) {
            anode::set_prefix(prefix);

            key_metas = new slot[cur_associativity * BUCKET_NUM]();
        }

        void traverse_for_pgm_resize(page_manager* old_pm, page_manager* new_pm,
                                     group_type resize_type) {
            size_t old_normal_pgid = normal_pgid;
            size_t old_special_pgid = special_pgid;
            size_t new_pgid = -1;

            if (resize_type == group_type::SPECIAL_GROUP) {
                new_pgid = new_pm->require_group_id(group_type::SPECIAL_GROUP);
                set_special_pgid(new_pgid);
            } else {
                new_pgid = new_pm->require_group_id(group_type::NORMAL_GROUP);
                set_normal_pgid(new_pgid);
            }

            page_manager_agent old_pm_agent = old_pm->get_page_manager_agent(old_normal_pgid, old_special_pgid);

            page_manager_agent new_pm_agent = 
                resize_type == group_type::SPECIAL_GROUP ? 
                new_pm->get_page_manager_agent(-1, new_pgid) : new_pm->get_page_manager_agent(new_pgid, -1);

            for (int i = 0; i != BUCKET_NUM; i++) {
                for (int j = 0; j != cur_associativity; j++) {
                    slot* s = get_slot(i, j);

                    // ignore the slot that not belong to current resize group
                    // type
                    if (s->is_empty() || get_group_type(s) != resize_type)
                        continue;

                    // get the content from old page_manager and write it to the
                    // new page_manager
                    s->set_slot(new_pm_agent.insert_element(old_pm_agent.get_content_pointer(s),
                                                        s->get_length(),
                                                        old_pm_agent.get_value(s)));
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

        inline int get_index(slot* s) { return s - key_metas; }

        inline int get_index(int bucketid, int slotid) {
            return bucketid * cur_associativity + slotid;
        }

        inline void set_normal_pgid (size_t new_normal_pgid) { normal_pgid = new_normal_pgid; }
        inline void set_special_pgid (size_t new_special_pgid) { special_pgid = new_special_pgid; }

        /* 
         * For eliminating the index update in expand_key_metas_space
         * we store the column-store-index in v2k instead of row-store-index
         */
        inline slot* get_column_store_slot(int column_store_index) {
            return key_metas +
                   (cur_associativity * (column_store_index % BUCKET_NUM)) +
                   column_store_index / BUCKET_NUM;
        }

        inline int get_column_store_index(slot* s) {
            return BUCKET_NUM * (get_index(s) % cur_associativity) +
                   (get_index(s) / cur_associativity);
        }

        /*---------function that changed key_metas layout---------*/
        /*------------------ 1. expand function------------------*/
        int expand_key_metas_space() {
            uint64_t sta = get_time();

            // Already max associativity
            // We cannot expand anymore, return -1
            if (cur_associativity == ASSOCIATIVITY) return -1;

            // Get the associativity we need, expand 2 times of cur_associativity
            unsigned int need_associativity = cur_associativity << 1;
            if (need_associativity > ASSOCIATIVITY) {
                need_associativity = ASSOCIATIVITY;
            }

            // Allocate a bigger memory for new key_metas
            slot* new_key_metas = new slot[need_associativity * BUCKET_NUM]();

            for (int i = 0; i != BUCKET_NUM; i++) {
                for (int j = 0; j != need_associativity; j++) {
                    slot* cur_new_slot =
                        new_key_metas + i * need_associativity + j;
                    if (j < cur_associativity) {
                        slot* cur_slot = key_metas + i * cur_associativity + j;
                        cur_new_slot->set_slot(*cur_slot);
                    } else {
                        cur_new_slot->set_slot(0, 0, 0, 0);
                    }
                }
            }

            // Switch the old key_metas to the new key_metas and release the old
            // key_metas
            delete[] key_metas;
            key_metas = new_key_metas;

            int ret_slotid = cur_associativity;
            // update current associativity
            cur_associativity = need_associativity;
            uint64_t end = get_time();

            expand_cost_time += end - sta;

            return ret_slotid;
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
        inline size_t get_another_bucketid(page_manager_agent& pm_agent, slot* s,
                                           size_t current_bucketid) {
            const char* key = pm_agent.get_content_pointer(s);
            size_t bucketid1 = hash(key, s->get_length(), 1) % BUCKET_NUM;
            size_t bucketid2 = hash(key, s->get_length(), 2) % BUCKET_NUM;
            return current_bucketid == bucketid1 ? bucketid2 : bucketid1;
        }

        // Return a empty slot_id in bucketid
        int cuckoo_hash(size_t bucketid, htrie_map<CharT, T>* hm) {
            cuckoohash_total_num++;
            uint64_t sta = get_time();

            // Set up the backup for recovery if the cuckoo hash fail
            slot* key_metas_backup = new slot[BUCKET_NUM * cur_associativity]();
            memcpy(key_metas_backup, key_metas,
                   BUCKET_NUM * cur_associativity * sizeof(slot));

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

            page_manager_agent pm_agent =
                hm->pm->get_page_manager_agent(normal_pgid, special_pgid);

            map<T, int> searchPoint_wait_2_be_update;
            for (int cuckoo_hash_time = 0; cuckoo_hash_time != BUCKET_NUM * ASSOCIATIVITY * DEFAULT_CUCKOO_HASH_RATIO;
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
                slot* cur_process_slot = get_slot(cur_process_bucketid, cur_process_slotid);

                /* Check that whether the cur_process_slot is anti-moved */
                // Get the another bucketid the cur_process_slot can be at
                int cur_kick_to_bucketid = get_another_bucketid(
                    pm_agent, cur_process_slot, cur_process_bucketid);
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
                searchPoint_wait_2_be_update[pm_agent.get_value(cur_process_slot)] =
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

                    cuckoohash_cost_time += get_time() - sta;

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
                   BUCKET_NUM * cur_associativity * sizeof(slot));

            delete[] key_metas_backup;
            delete extra_slot;

            cuckoohash_cost_time += get_time() - sta;

            // The cuckoo hash time exceeds Max_loop return -1 as slotid to
            // indicate cuckoo hash failed
            return -1;
        }

        /*----------------searching in hash_node----------------*/
        // return <found?, iterator>
        // iterator:    if slotid==-1, bucket is full
        //              if slotid!=-1, slotid is the insert position
        found_result find_in_bucket(size_t bucketid, const CharT* key,
                                    size_t keysize, page_manager_agent& pm_agent) {
            // find the hitted slot in hashnode
            for (int i = 0; i != cur_associativity; i++) {
                slot* target_slot = get_slot(bucketid, i);
                if (target_slot->is_empty())
                    return found_result(false, T(), bucketid, i);

                if (key_equal(key, keysize, pm_agent.get_content_pointer(target_slot),
                            target_slot->get_length()))
                    return found_result(true, pm_agent.get_value(target_slot), bucketid, i);
            }
            return found_result(false, T(), bucketid, -1);
        }

        found_result search_kv_in_hashnode(const CharT* key, size_t keysize,
                                       page_manager* pm) {
            page_manager_agent pm_agent = 
                    pm->get_page_manager_agent(normal_pgid, special_pgid);
            // if found the existed target in bucket1 or bucket2, just
            // return the iterator for being modified or read
            size_t bucket_id1 = hash(key, keysize, 1) % BUCKET_NUM;
            found_result res1 = find_in_bucket(bucket_id1, key, keysize, pm_agent);
            if (res1.is_founded()) {
                return res1;
            }

            size_t bucket_id2 = hash(key, keysize, 2) % BUCKET_NUM;
            found_result res2 = find_in_bucket(bucket_id2, key, keysize, pm_agent);
            if (res2.is_founded()) {
                return res2;
            }

            // if the code reach here it means the target doesn't exist
            // we try our best return the iterator with empty slot
            if (res1.is_bucket_full()) {
                return res2;
            } else
                return res1;
        }

        void insert_kv_in_hashnode(const CharT* key,
                                                 size_t keysize, htrie_map* hm,
                                                 T v, found_result fr) {
            size_t bucketid = fr.bucketid;
            int slotid = fr.slotid;
            page_manager_agent pm_agent = hm->pm->get_page_manager_agent(normal_pgid, special_pgid);

            if (slotid == -1 && 
                (slotid = expand_key_metas_space()) == -1 &&
                (slotid = cuckoo_hash(bucketid, hm)) == -1) {

                const string& prefix = this->anode::get_prefix();

                trie_node* new_parent =
                    hm->burst(burst_package(key_metas, BUCKET_NUM,
                                            cur_associativity, pm_agent),
                              this->anode::get_parent(), prefix);

                hm->access_kv_in_htrie_map(new_parent, key, keysize, v, false,
                                            prefix.data(), prefix.size());
                delete this;
                return;
            }

            // now the slotid cannot be -1 and slotid is lower than
            // ASSOCIATIVITY
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

            if(!pm_agent.try_insert(keysize)) {
                hm->pm->resize(get_group_type(keysize), hm);
                pm_agent = hm->pm->get_page_manager_agent(normal_pgid, special_pgid);
            }

            // call htrie-map function: write_kv_to_page ()
            // return a slot with position that element been written
            target_slot->set_slot(pm_agent.insert_element(key, keysize, v));

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
    class burst_package {
       private:
        page_manager_agent pm_agent_;
        vector<slot> elems_;

       public:
        burst_package(slot* elems, size_t bucket_num,
                      size_t associativity, page_manager_agent pm_agent)
            : pm_agent_(pm_agent) {
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

            page_manager_agent new_pm_agent = new_pm->get_page_manager_agent(n_group_id, s_group_id);
            for (auto& s : elems_) {
                // ignore the slot that not belong to current resize group
                // type
                if (get_group_type(&s) != resize_type) continue;

                // get the content from old page_manager and write it to the
                // new page_manager
                s.set_slot(new_pm_agent.insert_element(pm_agent_.get_content_pointer(&s),
                                    s.get_length(), pm_agent_.get_value(&s)));
            }
            pm_agent_.set_page_group(resize_type, new_pm_agent.get_page_group(resize_type));
        }

        slot operator[](int index) { return elems_[index]; }
        const slot top() const { return elems_.back(); }
        void pop() { elems_.pop_back(); }

        size_t size() {return elems_.size(); }
        page_manager_agent get_agent() const { return pm_agent_; }

        void print_bp() {
            cout << "-----------------\n";
            for(int i=0; i!= elems_.size(); i++){
                elems_[i].print_slot(pm_agent_);
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
                    char* key = pm_agent_.get_content_pointer(&(elems_[i]));
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
    trie_node* burst(burst_package bp, trie_node* orig_parent, const string &orig_prefix) {
        burst_total_counter++;

        pm->register_burst_package(&bp);

        // The return header
        trie_node* ret_trie_root = new trie_node(orig_parent, orig_prefix.data(), orig_prefix.size());

        if (orig_parent == nullptr)
            set_t_root(ret_trie_root);
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

                char* new_key = bp.get_agent().get_content_pointer(&s) + common_prefix_keysize;
                size_t length_left = s.get_length() - common_prefix_keysize;
                T v = bp.get_agent().get_value(&s);

                access_kv_in_htrie_map(parent, new_key, length_left, v, false, prefix.data(), prefix.size());

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
               private:
                friend class page_group;
                inline static unsigned int calc_align(unsigned int n, unsigned align) {
                    return ((n + align - 1) & (~(align - 1)));
                }
                unsigned int cur_pos;
                char* content;

               public:
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

                ~page() { if (content != nullptr) { free(content); } }
            };

            page* pages;
            int cur_page_id;
            bool is_special;

           public:
            page_group() : pages(nullptr), cur_page_id(-1), is_special(false) {}

            void init_pg(int page_number, bool spe) {
                is_special = spe;
                cur_page_id = 0;
                pages = new page[spe ? DEFAULT_SPECIAL_PAGE_NUMBER : DEFAULT_NORMAL_PAGE_NUMBER]();
                pages[0].init_page(spe ? DEFAULT_SPECIAL_PAGE_SIZE
                                       : DEFAULT_NORMAL_PAGE_SIZE);
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
                    (is_special ? DEFAULT_SPECIAL_PAGE_SIZE
                                : DEFAULT_NORMAL_PAGE_SIZE)) {
                    cur_page_id++;
                    pages[cur_page_id].init_page(
                        is_special ? DEFAULT_SPECIAL_PAGE_SIZE
                                   : DEFAULT_NORMAL_PAGE_SIZE);
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

            inline size_t get_cur_page_id(){
                return cur_page_id;
            }

            inline size_t get_max_page_id(){
                return is_special ? DEFAULT_SPECIAL_PAGE_NUMBER : DEFAULT_NORMAL_PAGE_NUMBER;
            }

            inline size_t get_max_per_page_size() {
                return is_special ? DEFAULT_SPECIAL_PAGE_SIZE : DEFAULT_NORMAL_PAGE_SIZE;
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
                delete []pages;
            }
        };
        
       private:
        page_group* normal_pg;
        page_group* special_pg;

        size_t n_size;
        size_t s_size;

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

        inline page_manager_agent get_page_manager_agent(int n_pg,
                                                         int s_pg) {
            return page_manager_agent(n_pg == -1 ? nullptr : normal_pg + n_pg,
                                        s_pg == -1 ? nullptr : special_pg + s_pg);
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
                if( *it == rm_bp_ptr){ 
                    notify_list.erase(it);
                    return;
                }
            }
        }

        void resize(group_type resize_type,
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

       private:
        void init_a_new_page_group(group_type init_type, size_t page_group_index) {
            if (init_type == group_type::SPECIAL_GROUP) {
                s_size++;
                special_pg[page_group_index].init_pg(DEFAULT_SPECIAL_PAGE_NUMBER, true);
                return;
            } else if (init_type == group_type::NORMAL_GROUP) {
                n_size++;
                normal_pg[page_group_index].init_pg(DEFAULT_NORMAL_PAGE_NUMBER, false);
                return;
            } else {
                cout << "undefined type!" << endl;
                assert(false);
                exit(0);
                return;
            }
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

        void notify_burst_package(page_manager *new_pm, group_type resize_type) {
            for (auto bp_ptr : notify_list) {
                bp_ptr->update_burst_package(new_pm, resize_type);
            }
        }
    };

    class SearchPoint {
       private:
        anode* node;
        int index;

       public:
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
                page_manager_agent pm_agent = pm->get_page_manager_agent(hnode->normal_pgid, hnode->special_pgid);
                res = res + string(pm_agent.get_content_pointer(sl),
                                   sl->get_length());
            }
            return res;
        }
    };

    class found_result {
       public:
        bool found;
        const T v;
        size_t bucketid;
        int slotid;

        found_result(bool f, T vv, size_t bid, int sid)
            : found(f), v(vv), bucketid(bid), slotid(sid) {}

        bool is_founded() { return found; }

        bool is_bucket_full() { return slotid == -1; }
    };

    uint32_t longest_string_size;
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
                        new hash_node(orig_tnode, new_prefix, pm);
                    orig_tnode->add_child(key[ref_pos - 1], current_node);
              }
            } 

          } break;
          case node_type::HASH_NODE: {
                hash_node* hnode = (hash_node*)current_node;
                found_result res = hnode->search_kv_in_hashnode(key + ref_pos,
                                                        key_size - ref_pos, pm);
                if (findMode) {
                    return std::pair<bool, T>(res.found, res.v);
                } else {
                    hnode->insert_kv_in_hashnode(key + ref_pos, key_size - ref_pos, this, v, res);
                    return std::pair<bool, T>(res.found, res.v);
                }
          } break;
          default:
            cout << "wrong type!";
            exit(0);
        }
      }

      // find a key in node's only value
      found_result res = current_node->search_kv_in_node();

      if (findMode) {
        return std::pair<bool, T>(res.found, res.v);
      } else {
        current_node->insert_value_in_node(string(prefix_key, prefix_key_size) + string(key, key_size), v, this);
        return std::pair<bool, T>(true, v);
      }
    }

    void set_t_root(anode* node) { t_root = node; }

    void set_searchPoint_index(T v, int index) { v2k[v].set_index(index); }

    void set_v2k(T v, anode* node, int index) {
        v2k[v] = SearchPoint(node, index);
    }

    // function for batch updating the searchPoints to v2k
    void apply_the_changed_searchPoint(map<T, int>& searchPoints) {
        for (auto it = searchPoints.begin(); it != searchPoints.end(); it++)
            set_searchPoint_index(it->first, it->second);
    }

    boost::unordered_map<T, SearchPoint> v2k;
    anode* t_root;
    page_manager *pm;

   public:
    htrie_map()
        : t_root(nullptr), pm(new page_manager()), longest_string_size(0) {
        std::cout << "SET UP GROWING-CUCKOOHASH-TRIE MAP\n";
        cout << "GROW_ASSOCIATIVITY\n";
        cout << "PM\n";

        t_root = new hash_node(nullptr, string(), pm, ASSOCIATIVITY);
    }

    // TODO deconstructor

    /*---------------external accessing interface-------------------*/
    // TODO: adapt to wukong's interface
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
        pm->resize(group_type ::NORMAL_GROUP, this, 0);
        pm->resize(group_type ::SPECIAL_GROUP, this, 0);
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