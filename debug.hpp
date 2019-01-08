#include <unistd.h>
#include <fstream>
#include <iostream>
#include <vector>

namespace myTrie {
namespace debuging {
using std::cout;
using std::endl;

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
double last_time_vm_usage = 0;
double last_time_resident_set = 0;

void clear_process_mem_usage() {
    using std::ifstream;
    using std::ios_base;
    using std::string;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
        tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
        stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
        starttime >> vsize >> rss;  // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                        1024;  // in case x86-64 is configured to use 2MB pages

    unsigned long now_vm_usage = vsize;
    long now_resident_set = rss * page_size_kb;

    last_time_vm_usage = now_vm_usage;
    last_time_resident_set = now_resident_set;
}

void process_mem_usage(double& vm_usage, double& resident_set) {
    using std::ifstream;
    using std::ios_base;
    using std::string;

    vm_usage = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
        tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
        stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
        starttime >> vsize >> rss;  // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                        1024;  // in case x86-64 is configured to use 2MB pages
    vm_usage = (vsize - last_time_vm_usage) / 1024.0;
    resident_set = rss * page_size_kb - last_time_resident_set;
    std::cout << "cur proc vm_usage: " << vm_usage << " res: " << resident_set
              << std::endl;
}
}  // namespace debuging
}  // namespace myTrie

namespace myTrie {
namespace debuging {
using std::cout;
using std::endl;
uint64_t h_n;
uint64_t t_n;
uint64_t m_n;

uint64_t hash_node_mem;
uint64_t trie_node_mem;
uint64_t multi_node_mem;
uint64_t v2k_mem;

uint64_t child_first_char_mem;

size_t hashnode_load = 0;
#ifdef TEST_CUCKOOHASH
size_t hashnode_max_load = Max_slot_num / 2;
size_t hashnode_min_load = Max_slot_num / 2;
#endif

#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
size_t hashnode_max_load = Max_slot_num / 2;
size_t hashnode_min_load = Max_slot_num / 2;
size_t hashnode_total_slot_num = 0;
#endif

size_t total_pass_trie_node_num = 0;
size_t myCount = 0;

template <typename CharT, typename T>
void print_tree_construct_v2k(
    map<T, class myTrie::htrie_map<CharT, T>::SearchPoint>& v2k) {
    size_t counter_sz = v2k.size();
    size_t entry_sz =
        sizeof(T) + sizeof(class myTrie::htrie_map<CharT, T>::SearchPoint);
    cout << "sizeof(class myTrie::htrie_map<CharT, T>::SearchPoint): "
         << sizeof(class myTrie::htrie_map<CharT, T>::SearchPoint) << endl;

    cout << "sizeof(int32_t): " << sizeof(int32_t) << endl;
    cout << "sizeof(void*): " << sizeof(void*) << endl;

    v2k_mem = counter_sz * entry_sz;
}

template <typename CharT, typename T>
void print_tree_construct(class myTrie::htrie_map<CharT, T>::anode* root,
                          size_t depth = 0) {
    using my = myTrie::htrie_map<CharT, T>;
    if (root == nullptr) {
        return;
    }
    // print bucket
    if (root->is_hash_node()) {
        h_n++;

        class my::hash_node* cur_hash_node = (class my::hash_node*)root;
#ifdef TEST_CUCKOOHASH
        hash_node_mem += sizeof(cur_hash_node);

        // basic bucket cost
        uint32_t bucket_mem = sizeof(class my::hash_node::slot) * Max_slot_num;
        hash_node_mem += bucket_mem;

        // page mem cost
        size_t pages_cost = Max_bytes_per_kv * cur_hash_node->pages.size();
        hash_node_mem += pages_cost;

        h_n++;

        size_t current_elem_num = cur_hash_node->elem_num;
        hashnode_load += current_elem_num;
        if (current_elem_num > hashnode_max_load) {
            hashnode_max_load = current_elem_num;
        }

        if (current_elem_num < hashnode_min_load) {
            hashnode_min_load = current_elem_num;
        }

        total_pass_trie_node_num += cur_hash_node->elem_num * depth;

        if (cur_hash_node->haveValue) {
            total_pass_trie_node_num += depth;
        }
#endif
#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
        hash_node_mem += sizeof(cur_hash_node);
        size_t entry_sz = 0;
        size_t counter_sz = 0;
        // skip the data structure that used to burst
        // #ifdef IMPROVE_BURST
        //         counter_sz = (cur_hash_node->element_num_of_1st_char).size();
        //         entry_sz = sizeof(CharT) + sizeof(uint16_t);
        // #endif
        // #ifdef BURST_IN_ADVANCE
        //         counter_sz = (cur_hash_node->element_num_of_1st_char).size();
        //         entry_sz = sizeof(CharT) + sizeof(class
        //         my::hash_node::triple);
        // #endif

        hash_node_mem += counter_sz * entry_sz;
        child_first_char_mem += counter_sz * entry_sz;

        // basic bucket cost
#ifdef GROW_BUCKET
        uint32_t bucket_mem = sizeof(class my::hash_node::slot) *
                              (cur_hash_node->cur_bucket) * Associativity;
#else
        uint32_t bucket_mem = sizeof(class my::hash_node::slot) *
                              (cur_hash_node->cur_associativity) * Bucket_num;
#endif

        hash_node_mem += bucket_mem;

        // page mem cost
        size_t pages_cost = Max_bytes_per_kv * cur_hash_node->pages.size();
        hash_node_mem += pages_cost;

        size_t current_elem_num = cur_hash_node->elem_num;
        hashnode_load += current_elem_num;
#ifdef GROW_BUCKET
        hashnode_total_slot_num += (cur_hash_node->cur_bucket) * Associativity;
#else
        hashnode_total_slot_num +=
            (cur_hash_node->cur_associativity) * Bucket_num;
#endif

        if (current_elem_num > hashnode_max_load) {
            hashnode_max_load = current_elem_num;
        }

        if (current_elem_num < hashnode_min_load) {
            hashnode_min_load = current_elem_num;
        }

        total_pass_trie_node_num += cur_hash_node->elem_num * depth;
#endif
#ifdef TEST_HAT
        hash_node_mem += sizeof(cur_hash_node);

        std::vector<class my::array_bucket> buckets = cur_hash_node->kvs;
        // basic bucket cost
        uint32_t bucket_num = buckets.size();
        uint32_t bucket_mem = sizeof(class my::array_bucket) * bucket_num;
        hash_node_mem += bucket_mem;

        uint32_t bucket_buffer_mem = 0;
        // bucket buffer cost
        for (auto it = buckets.begin(); it != buckets.end(); it++) {
            bucket_buffer_mem += it->buffer_size;
        }
        hash_node_mem += bucket_buffer_mem;

        hashnode_load += cur_hash_node->element_count * depth;

        total_pass_trie_node_num += cur_hash_node->element_count * depth;

        if (cur_hash_node->haveValue) total_pass_trie_node_num += depth;
        h_n++;
#endif

    } else if (root->is_trie_node()) {
        t_n++;
        class myTrie::htrie_map<CharT, T>::trie_node* cur_trie_node =
            (class myTrie::htrie_map<CharT, T>::trie_node*)root;
        trie_node_mem += sizeof(cur_trie_node);
#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
        if (cur_trie_node->prefix_len != 0) {
            trie_node_mem += cur_trie_node->prefix_len;
        }
#endif

#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
        std::map<CharT, class myTrie::htrie_map<CharT, T>::trie_node*> childs =
            cur_trie_node->trie_node_childs;
#else
        std::map<CharT, class myTrie::htrie_map<CharT, T>::trie_node*> childs =
            cur_trie_node->childs;
#endif

        if (childs.size() == 0) {
#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
            print_tree_construct<CharT, T>(cur_trie_node->hash_node_child,
                                           depth + 1);
#else
            print_tree_construct<CharT, T>(cur_trie_node->onlyHashNode,
                                           depth + 1);
#endif
        } else {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_construct<CharT, T>(it->second, depth + 1);
            }
        }
    } else {
#ifdef SHRINK_TEST_GROWCUCKOOHASH
        if (root->is_multi_node()) {
            m_n++;

            class myTrie::htrie_map<CharT, T>::multi_node* cur_multi_node =
                (class myTrie::htrie_map<CharT, T>::multi_node*)root;
            multi_node_mem += sizeof(cur_multi_node);
            std::map<string, class myTrie::htrie_map<CharT, T>::anode*> childs =
                cur_multi_node->childs_;
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_construct<CharT, T>(it->second, depth + 1);
            }
            return;
        }
#endif
        std::cout << "node is not trie nor hash node\n";
        exit(0);
    }
    return;
}

void clear_num() {
    t_n = 0;
    h_n = 0;
    hash_node_mem = 0;
    trie_node_mem = 0;
    multi_node_mem = 0;
    v2k_mem = 0;
    child_first_char_mem = 0;

    hashnode_load = 0;
#ifdef TEST_CUCKOOHASH
    hashnode_max_load = Max_slot_num / 2;
    hashnode_min_load = Max_slot_num / 2;
#endif
#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
    hashnode_max_load = Max_slot_num / 2;
    hashnode_min_load = Max_slot_num / 2;
    hashnode_total_slot_num = 0;
#endif

    total_pass_trie_node_num = 0;
    return;
}

template <typename CharT, typename T>
double print_res() {
#ifdef TEST_CUCKOOHASH
    std::cout << "WE ARE TESTING CUCKOOHASH_IMPL NOW\n";
#endif
#ifdef TEST_GROWCUCKOOHASH
    std::cout << "WE ARE TESTING GROWING_CUCKOOHASH_IMPL NOW\n";
#endif
#ifdef SHRINK_TEST_GROWCUCKOOHASH
    std::cout << "WE ARE TESTING SHRINKING GROWING_CUCKOOHASH_IMPL NOW\n";
#endif
#ifdef TEST_HAT
    std::cout << "WE ARE TESTING HAT-TRIE NOW\n";
#endif

    std::cout << "trie_node: " << t_n << " hash_node: " << h_n
              << " multi_node: " << m_n << std::endl;
    std::cout << "trie_node_mem: " << trie_node_mem
              << " hash_node_mem: " << hash_node_mem
              << " multi_node_mem: " << multi_node_mem
              << " child_first_char_mem: " << child_first_char_mem
              << " v2k_mem: " << v2k_mem << std::endl;

    using my = typename myTrie::htrie_map<CharT, T>;

    std::cout << "total: ,"
              << (hash_node_mem + trie_node_mem + multi_node_mem + v2k_mem +
                  child_first_char_mem) /
                     1024 / 1024
              << ", mb" << std::endl;
    double ret_v;
    ret_v = (hash_node_mem + trie_node_mem + multi_node_mem + v2k_mem +
             child_first_char_mem) /
            1024 / 1024;

    // todo
    // std::cout << "trienode size:" << sizeof(typename my::trie_node)
    //           << std::endl;
    // std::cout << "hashnode size:" << sizeof(typename my::hash_node)
    //           << std::endl;
    // std::cout << "arraybucket size:" << sizeof(typename my::array_bucket)
    //           << std::endl;
    // std::cout << "htrie_map size:"
    //           << sizeof(typename myTrie::htrie_map<CharT, T>) << std::endl;
    // std::cout
    //     << "htrie_map::v2k size:"
    //     << sizeof(
    //            std::map<T, typename myTrie::htrie_map<CharT,
    //            T>::SearchPoint>)
    //     << std::endl;
    // std::cout << "arraybucket::kvs size:"
    //           << sizeof(std::vector<typename my::array_bucket>) << std::endl;

    return ret_v;
}
}  // namespace debuging
}  // namespace myTrie

uint64_t get_usec() {
    struct timespec tp;
    /* POSIX.1-2008: Applications should use the clock_gettime() function
       instead of the obsolescent gettimeofday() function. */
    /* NOTE: The clock_gettime() function is only available on Linux.
       The mach_absolute_time() function is an alternative on OSX. */
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return ((tp.tv_sec * 1000 * 1000) + (tp.tv_nsec / 1000));
}

#include <sys/sysinfo.h>
uint64_t getLftMem() {
    struct sysinfo inf;
    sysinfo(&inf);

    return inf.mem_unit *
           (inf.totalram + inf.totalswap - inf.freeram - inf.freeswap);
}