#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "test_which.hpp"
// #include "unified_impl/5_prototype_cuckoo_grow_buc_shrink.hpp"
// #include "unified_impl/4_prototype_cuckoo_grow_ass_shrink.hpp"
// #include "unified_impl/3_prototype_cuckoo_shrink.hpp"
// #include "unified_impl/2_prototype.hpp"
// #include "unified_impl/1_tessil_hat_impl.hpp"

#include "debug.hpp"

#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>

#include <sys/sysinfo.h>

using namespace std;
int main() {
    myTrie::htrie_map<char, uint32_t> hm(4, 31);
    map<string, uint32_t> m1;
    map<uint32_t, string> m2;

#ifdef TEST_YAGO
    string testing_dataset = "dataset/id_yago/cut_str_normal";
#else
    string testing_dataset = "dataset/str_normal";
#endif

    cout << "testing file: " << testing_dataset << endl;

    fstream f(testing_dataset);
    string url;
    uint32_t v;
    uint64_t sta = get_time();
    while (f >> url >> v) {
        hm.insertKV(url, v);
        m1[url] = v;
        m2[v] = url;
    }
    uint64_t end = get_time();
    std::cout << "loaded 3 map in " << (end - sta) / 1000 << " ms" << std::endl;

    // myTrie::debuging::print_tree_construct<char, uint32_t>(hm.t_root);
    // double mem_cal_inside = myTrie::debuging::print_res<char, uint32_t>();

#ifdef SHRINK_TEST_GROWCUCKOOHASH
    hm.shrink();
#endif

#ifdef TEST_GROWCUCKOOHASH
    cout << "expand cost time: " << expand_cost_time << endl;
#endif

    uint64_t count = m1.size();
    uint64_t searchByK_count_good = 0;
    uint64_t searchByK_count_bad = 0;
    uint64_t searchByV_count_good = 0;
    uint64_t searchByV_count_bad = 0;
    vector<uint32_t> vv;

    for (auto it = m1.begin(); it != m1.end(); it++) {
        if (hm.searchByKey(it->first) == it->second) {
            // std::cout << "good\n";
            // std::cout << "ans: " << it->second << std::endl;
            // std::cout << "got " << hm.searchByKey(it->first).second
            //           << std::endl;
            searchByK_count_good++;
        } else {
            std::cout << "wrong\n";
            std::cout << "ans: " << it->second << std::endl;
            vv.push_back(it->second);
            std::cout << "got " << hm.searchByKey(it->first) << std::endl;
            searchByK_count_bad++;
        }
    }
    std::cout << "-------------------\n";

    for (auto it = m2.begin(); it != m2.end(); it++) {
        if (hm.searchByValue(it->first) == it->second) {
            // std::cout << "good\n";
            // std::cout << "ans: " << it->second << std::endl;
            // std::cout << "got " << hm.searchByValue(it->first) << std::endl;
            searchByV_count_good++;
        } else {
            std::cout << "wrong\n";
            std::cout << "ans: " << it->second << std::endl;
            std::cout << "got " << hm.searchByValue(it->first) << std::endl;
            searchByV_count_bad++;
        }
    }
    std::cout << "-------------------\n";

    cout << "total: " << count << endl;
    cout << "check key: good:" << searchByK_count_good
         << " bad:" << searchByK_count_bad << std::endl;
    cout << "check value: good:" << searchByV_count_good
         << " bad:" << searchByV_count_bad << std::endl;
    cout << "finish\n";
}