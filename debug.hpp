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

#include "myTrie.hpp"

namespace myTrie {
namespace debuging {
using std::cout;
using std::endl;
uint64_t h_n;
uint64_t t_n;
uint64_t hash_node_mem;
uint64_t trie_node_mem;

size_t hashnode_load = 0;
size_t hashnode_max_load = Max_slot_num / 2;
size_t hashnode_min_load = Max_slot_num / 2;

size_t total_pass_trie_node_num = 0;
size_t myCount = 0;

template <typename CharT, typename T>
void print_tree_construct(typename myTrie::htrie_map<CharT, T>::anode* root,
                          size_t depth = 0) {
    using my = typename myTrie::htrie_map<CharT, T>;
    if (root == nullptr) {
        return;
    }
    // print bucket
    if (root->isHashNode()) {
        hash_node_mem += sizeof(root);

        // basic bucket cost
        uint32_t bucket_mem =
            sizeof(typename my::hash_node::slot) * Max_slot_num;
        hash_node_mem += bucket_mem;

        // page mem cost
        size_t pages_cost =
            Max_bytes_per_kv * ((typename my::hash_node*)root)->pages.size();
        hash_node_mem += pages_cost;

        h_n++;

        size_t current_elem_num = ((typename my::hash_node*)root)->elem_num;
        hashnode_load += current_elem_num;
        if (current_elem_num > hashnode_max_load) {
            hashnode_max_load = current_elem_num;
        }

        if (current_elem_num < hashnode_min_load) {
            hashnode_min_load = current_elem_num;
        }

        total_pass_trie_node_num +=
            ((typename my::hash_node*)root)->elem_num * depth;

        if (((typename my::hash_node*)root)->haveValue) {
            total_pass_trie_node_num += depth;
        }

    } else if (root->isTrieNode()) {
        trie_node_mem += sizeof(root);

        std::map<CharT, typename myTrie::htrie_map<CharT, T>::trie_node*>
            childs = ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                         ->childs;
        if (childs.size() == 0) {
            print_tree_construct<CharT, T>(
                ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                    ->getOnlyHashNode(),
                depth + 1);
        } else {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                t_n++;
                print_tree_construct<CharT, T>(it->second, depth + 1);
            }
        }
    } else {
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
    hashnode_load = 0;
    hashnode_max_load = Max_slot_num / 2;
    hashnode_min_load = Max_slot_num / 2;
    total_pass_trie_node_num = 0;
    return;
}

template <typename CharT, typename T>
double print_res() {
    std::cout << "trie_node: " << t_n << " hash_node: " << h_n << std::endl;
    std::cout << "trie_node_mem: " << trie_node_mem
              << " hash_node_mem: " << hash_node_mem << std::endl;
    using my = typename myTrie::htrie_map<CharT, T>;

    std::cout << "total: ," << (hash_node_mem + trie_node_mem) / 1024 / 1024
              << ", mb" << std::endl;
    double ret_v;
    ret_v = (hash_node_mem + trie_node_mem) / 1024 / 1024;

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