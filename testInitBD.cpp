#include <fstream>
#include <iostream>
#include <map>
#include <string>

// #include "impl_growing_cuckoo_hash_shrink/grow_cuckoo_trie_ass_shrnk_improve_burst_inadvance.hpp"
#include "impl_growing_cuckoo_hash_shrink/grow_cuckoo_trie_ass_shrnk_improve_burst.hpp"

#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>

#include <sys/sysinfo.h>

static uint64_t get_usec() {
    struct timespec tp;
    /* POSIX.1-2008: Applications should use the clock_gettime() function
       instead of the obsolescent gettimeofday() function. */
    /* NOTE: The clock_gettime() function is only available on Linux.
       The mach_absolute_time() function is an alternative on OSX. */
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return ((tp.tv_sec * 1000 * 1000) + (tp.tv_nsec / 1000));
}

using namespace std;
int main() {
    myTrie::htrie_map<char, uint32_t> hm(4, 5);
    // string testing_dataset = "dataset/str_normal";
    string testing_dataset = "dataset/id_yago/cut_str_normal";
    fstream f(testing_dataset);

    string url;
    uint32_t v;

    uint64_t total_cost_time = 0;
    while (f >> url >> v) {
        uint64_t sta = get_usec();
        hm.insertKV(url, v);
        uint64_t end = get_usec();
        total_cost_time += end - sta;
    }
    f.close();

    map<uint32_t, uint32_t> counting_reverse;
    for (auto it = insert_time_counting.begin();
         it != insert_time_counting.end(); it++) {
        counting_reverse[it->second]++;
    }
    for (auto it = counting_reverse.begin(); it != counting_reverse.end();
         it++) {
        cout << "insert time: " << it->first << " number: " << it->second
             << endl;
    }

    cout << "total cost time: " << total_cost_time << endl;
    cout << "expand cost time: " << expand_cost_time << endl;
    cout << "rehash_cost_time cost time: " << rehash_cost_time << endl;
    cout << "burst_cost_time cost time: " << burst_cost_time << endl;
    cout << "-burst_preprocess_time cost time: " << burst_preprocess_time
         << endl;
    cout << "insert_in_trienode_cost_time cost time: "
         << insert_in_trienode_cost_time << endl;
    cout << "insert_in_hashnode_cost_time cost time: "
         << insert_in_hashnode_cost_time << endl;
    cout << "non_shrink_find_hashnode_in_trienode_cost_time cost time: "
         << non_shrink_find_hashnode_in_trienode_cost_time << endl;
    cout << "non_shrink_find_hashnode_in_hashnode_cost_time cost time: "
         << non_shrink_find_hashnode_in_hashnode_cost_time << endl;
    cout << "time_is_zero_time: " << time_is_zero_time << endl;

    cout << "recal_element_num_of_1st_char_counter: " << recal_element_num_of_1st_char_counter << endl;
    cout << "burst_total_counter: " << burst_total_counter << endl;

    f.open(testing_dataset);
    map<string, uint32_t> helper_map;
    while (f >> url >> v) {
        helper_map[url] = v;
    }
    cout << "helper_map.size(): " << helper_map.size() << endl;

    non_shrink_find_hashnode_in_trienode_cost_time = 0;
    non_shrink_find_hashnode_in_hashnode_cost_time = 0;

    for (auto it = helper_map.begin(); it != helper_map.end(); it++) {
        if (hm.searchByKey(it->first) == it->second) {
        } else {
            cout << "wrong!";
        }
    }

    cout << "non_shrink_find_hashnode_in_trienode_cost_time cost time: "
         << non_shrink_find_hashnode_in_trienode_cost_time << endl;
    cout << "non_shrink_find_hashnode_in_hashnode_cost_time cost time: "
         << non_shrink_find_hashnode_in_hashnode_cost_time << endl;

    uint64_t sta1 = get_usec();
    hm.shrink();
    uint64_t end1 = get_usec();

    for (auto it = helper_map.begin(); it != helper_map.end(); it++) {
        if (hm.searchByKey(it->first) == it->second) {
        } else {
            cout << "wrong!";
        }
    }

    cout << "shrink_total cost time: " << end1 - sta1 << endl;
    cout << "shrink_find_hashnode_in_multinode_cost_time cost time: "
         << shrink_find_hashnode_in_multinode_cost_time << endl;
    cout << "shrink_find_hashnode_in_trienode_cost_time cost time: "
         << shrink_find_hashnode_in_trienode_cost_time << endl;
    cout << "shrink_find_hashnode_in_hashnode_cost_time cost time: "
         << shrink_find_hashnode_in_hashnode_cost_time << endl;
}