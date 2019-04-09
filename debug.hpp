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

void process_mem_usage(std::string part_name, double& vm_usage, double& resident_set) {
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
    std::cout << part_name << " : vm_usage: " << vm_usage << " res: " << resident_set
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

uint64_t child_representation_mem;

size_t hashnode_max_load = Max_slot_num / 2;
size_t hashnode_min_load = Max_slot_num / 2;

size_t hashnode_total_slot_num = 0;
size_t hashnode_total_element = 0;

size_t total_pass_trie_node_num = 0;
size_t myCount = 0;

size_t byte_used_in_page = 0;
size_t byte_pages_have = 0;
size_t page_metadata_mem = 0;
size_t total_normal_page_number = 0;
size_t total_special_page_number = 0;

size_t prefix_extra_mem = 0;

size_t hashnode_keymetas_mem = 0;

/*
traverse for
    Memory: trie_node_mem, hash_node_mem, bucket_mem, hashnode_keymetas_mem, prefix_extra_mem
    Number: hashnode_total_slot_num, hashnode_total_element, total_pass_trie_node_num
    t_n, h_n
 */
template <typename CharT, typename T>
void traverse_trie(class myTrie::htrie_map<CharT, T>::anode* root,
                          size_t depth = 0) {
    using my = myTrie::htrie_map<CharT, T>;
    if (root == nullptr) {
        return;
    }
    // print bucket
    if (root->is_hash_node()) {
        h_n++;

        class my::hash_node* cur_hash_node = (class my::hash_node*)root;
        hash_node_mem += sizeof(cur_hash_node);

/*---------------------- bucket mem cost ----------------------*/
#ifndef TEST_HAT

#ifdef TEST_GROW_CUCKOOHASH_BUC
        uint32_t bucket_mem = sizeof(class my::hash_node::slot) *
                              (cur_hash_node->cur_bucket) * Associativity;
        hashnode_total_slot_num += (cur_hash_node->cur_bucket) * Associativity;

#endif
#ifdef TEST_GROW_CUCKOOHASH_ASS
        uint32_t bucket_mem = sizeof(class my::slot) *
                              (cur_hash_node->cur_associativity) * Bucket_num;
        hashnode_total_slot_num +=
            (cur_hash_node->cur_associativity) * Bucket_num;

#endif
#if (defined TEST_CUCKOOHASH) || (defined TEST_HASH)
        uint32_t bucket_mem =
            sizeof(class my::hash_node::slot) * Associativity * Bucket_num;
        hashnode_total_slot_num += Associativity * Bucket_num;

#endif
        // add the bucket mem to hash_node_mem
        hash_node_mem += bucket_mem;
#else
        // add the bucket mem to hash_node_mem
        std::vector<class my::hash_node::array_bucket>& buckets =
            cur_hash_node->kvs;

        for (int i = 0; i != buckets.size(); i++) {
            size_t current_bucket_mem =
                sizeof(class my::hash_node::array_bucket);

            hash_node_mem += buckets[i].buffer_size;
            hash_node_mem += sizeof(class my::hash_node::array_bucket);
        }

        hashnode_total_slot_num += Associativity * Bucket_num;
#endif

#if (defined TEST_HASH) || (defined TEST_CUCKOOHASH) || \
    (defined TEST_GROW_CUCKOOHASH_ASS) || (defined TEST_GROW_CUCKOOHASH_BUC)
        hashnode_keymetas_mem += bucket_mem;
#endif

        // load ratio calculation
        size_t current_elem_num = cur_hash_node->elem_num;

        hashnode_total_element += current_elem_num;

        if (current_elem_num > hashnode_max_load) {
            hashnode_max_load = current_elem_num;
        }

        if (current_elem_num < hashnode_min_load) {
            hashnode_min_load = current_elem_num;
        }

        // total parent calculation
        total_pass_trie_node_num += cur_hash_node->elem_num * depth;
    } else if (root->is_trie_node()) {
        t_n++;
        class myTrie::htrie_map<CharT, T>::trie_node* cur_trie_node =
            (class myTrie::htrie_map<CharT, T>::trie_node*)root;
        trie_node_mem += sizeof(cur_trie_node);

        if (cur_trie_node->prefix_len != 0) {
            trie_node_mem += cur_trie_node->prefix_len;
        }

        vector<pair<CharT, class myTrie::htrie_map<CharT, T>::trie_node*>>
            childs;

         cur_trie_node->trie_node_childs.get_childs_with_char(childs);

        prefix_extra_mem +=
            sizeof(CharT*) + sizeof(uint16_t) + cur_trie_node->prefix_len;

        if (childs.size() == 0) {
            traverse_trie<CharT, T>(cur_trie_node->get_hash_node_child(),
                                           depth + 1);
        } else {
            
            for (auto it = childs.begin(); it != childs.end(); it++) {
                traverse_trie<CharT, T>(it->second, depth + 1);
            }
        }

         size_t cur_child_representation_mem =
             cur_trie_node->trie_node_childs.get_childs_representation_mem();
         child_representation_mem += cur_child_representation_mem;
         trie_node_mem += cur_child_representation_mem;

    } else {
        if (root->is_multi_node()) {
            m_n++;

            class myTrie::htrie_map<CharT, T>::multi_node* cur_multi_node =
                (class myTrie::htrie_map<CharT, T>::multi_node*)root;
            multi_node_mem += sizeof(cur_multi_node);
            std::map<string, class myTrie::htrie_map<CharT, T>::anode*> childs =
                cur_multi_node->childs_;
            for (auto it = childs.begin(); it != childs.end(); it++) {
                traverse_trie<CharT, T>(it->second, depth + 1);
            }
            return;
        }
        std::cout << "node is not trie nor hash node\n";
        exit(0);
    }
    return;
}


template <typename CharT, typename T>
void scan_tree(class myTrie::htrie_map<CharT, T> &hm){
    // traverse trie to get memory,number information
    traverse_trie<CharT, T>(hm.t_root);

    // scan normal
    for (int i = 0; i != hm.normal_pages.size(); i++) {
        size_t page_used = hm.normal_pages[i].cur_pos;
        size_t page_have = Max_bytes_per_kv;
        
        byte_pages_have += page_have;
        byte_used_in_page += page_used;

        total_normal_page_number++;
        page_metadata_mem += sizeof(class myTrie::htrie_map<CharT, T>::page);
    }

    // scan special
    for (int i = 0; i != hm.special_pages.size(); i++) {
        size_t page_used = hm.special_pages[i].cur_pos;
        size_t page_have = DEFAULT_SPECIAL_Max_bytes_per_kv;
        
        byte_pages_have += page_have;
        byte_used_in_page += page_used;

        total_special_page_number++;
        page_metadata_mem += sizeof(class myTrie::htrie_map<CharT, T>::page);
    }

    // value to searchPoint scanning
    size_t counter_sz = hm.v2k.size();
    size_t entry_sz =
        sizeof(T) + sizeof(class myTrie::htrie_map<CharT, T>::SearchPoint);
    v2k_mem = counter_sz * entry_sz;

}

void clear_num() {
    t_n = 0;
    h_n = 0;
    hash_node_mem = 0;
    trie_node_mem = 0;
    multi_node_mem = 0;
    v2k_mem = 0;

    child_representation_mem = 0;

    hashnode_total_element = 0;

    hashnode_max_load = Max_slot_num / 2;
    hashnode_min_load = Max_slot_num / 2;
    hashnode_total_slot_num = 0;

    total_pass_trie_node_num = 0;

    byte_used_in_page = 0;
    byte_pages_have = 0;
    page_metadata_mem = 0;
    total_normal_page_number = 0;
    total_special_page_number = 0;

    prefix_extra_mem = 0;
    hashnode_keymetas_mem = 0;
    return;
}

template <typename CharT, typename T>
double print_res() {
    cout << "---------debug result printing-------------" << endl;
#ifdef TEST_HAT
    cout << "unified_impl/1_tessil_hat_impl.hpp" << endl;
#endif
#ifdef TEST_HASH
    cout << "unified_impl/2_prototype_shrink.hpp" << endl;
#endif
#ifdef TEST_CUCKOOHASH
    cout << "unified_impl/3_prototype_cuckoo_shrink.hpp" << endl;
#endif
#ifdef TEST_GROW_CUCKOOHASH_ASS
    cout << "unified_impl/4_prototype_cuckoo_grow_ass_shrink.hpp" << endl;
#endif
#ifdef TEST_GROW_CUCKOOHASH_BUC
    cout << "unified_impl/5_prototype_cuckoo_grow_buc_shrink.hpp" << endl;
#endif

    cout << "trie_node: " << t_n << " hash_node: " << h_n
         << " multi_node: " << m_n << std::endl;
    cout << "trie_node_mem: " << trie_node_mem
         << " hash_node_mem: " << hash_node_mem
         << " multi_node_mem: " << multi_node_mem
         << " v2k_mem: " << v2k_mem 
         << " byte_pages_have: " << byte_pages_have << std::endl;
    cout << "page load: "
         << (double)byte_used_in_page / (double)byte_pages_have * 100 << "%"
         << endl;

    using my = typename myTrie::htrie_map<CharT, T>;

    cout << "total: ,"
         << (hash_node_mem + trie_node_mem + multi_node_mem + v2k_mem) / 1024 /
                1024
         << ", mb" << std::endl;
    double ret_v;
    ret_v = (hash_node_mem + trie_node_mem + multi_node_mem + v2k_mem) / 1024 /
            1024;
    cout << "prefix_extra_mem: " << prefix_extra_mem << endl;
    cout << "Kmetas avg load: "
         << (double)hashnode_total_element / (double)hashnode_total_slot_num *
                100
         << "%" << endl;

    cout << "\n-----BREAKDOWN-----" << endl;
    cout << "byte_pages_have: " << byte_pages_have << endl;
    cout << "hashnode_keymetas_mem: " << hashnode_keymetas_mem << endl;
    cout << "child_representation_mem: " << child_representation_mem << endl;
    cout << "page_metadata_mem: " << page_metadata_mem << endl;
    
    return ret_v;
}
}  // namespace debuging
}  // namespace myTrie

// #include <sys/sysinfo.h>
// uint64_t getLftMem() {
//     struct sysinfo inf;
//     sysinfo(&inf);

//     return inf.mem_unit *
//            (inf.totalram + inf.totalswap - inf.freeram - inf.freeswap);
// }