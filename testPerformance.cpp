bool test_and_print_wrong_test = true;
bool manually_test = false;

#include "test_which.hpp"

// #include "myTrie.hpp"
// debug
#include <fstream>
#include "debug.hpp"

#include <boost/unordered_map.hpp>
#include <map>
#include <vector>
using std::fstream;

#include <stdint.h>
#include <sys/time.h>

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
    std::cout << "starto!\n";

    std::vector<std::pair<size_t, size_t>> configs;
    vector<size_t> elem_per_buck;
    vector<size_t> bucket_nums;

    // analyse by elem_per_bucket
    // associativity should be 4, 8 for alignment
    elem_per_buck.push_back(4);
    elem_per_buck.push_back(8);

    // analyse by bucket_num
    bucket_nums.push_back(5);
    bucket_nums.push_back(11);
    bucket_nums.push_back(31);
    bucket_nums.push_back(59);
    bucket_nums.push_back(71);
    bucket_nums.push_back(97);
    bucket_nums.push_back(101);

#ifdef TEST_YAGO
    string testing_dataset = "dataset/id_yago/cut_str_normal";
#else
    string testing_dataset = "dataset/str_normal";
#endif

    cout << "testing file: " << testing_dataset << endl;

    // configs: pair<bucket_num, elem_per_bucket>
    // stable bucket_num
    for (auto itt = elem_per_buck.begin(); itt != elem_per_buck.end(); itt++) {
        for (auto it = bucket_nums.begin(); it != bucket_nums.end(); it++) {
            configs.push_back(std::pair<size_t, size_t>(*itt, *it));
        }
    }

    fstream ff1("result", std::ios::out | std::ios::app);

    //------------testing--------------
    std::string url;
    uint32_t v;
    uint64_t staUm;
    uint64_t endUm;
    uint64_t staTm;
    uint64_t endTm;

    std::fstream f(testing_dataset);
    boost::unordered_map<string, uint32_t> m1;
    boost::unordered_map<uint32_t, string> m2;

    staUm = get_usec();
    uint64_t startUsedMemUm = getLftMem();

    uint32_t count = 0;
    while (f >> url >> v) {
        m1[url] = v;
        m2[v] = url;
        count++;
    }
    uint64_t endUsedMemUm = getLftMem();
    endUm = get_usec();
    f.close();
    cout << "unordered_map use time: usec: \t\t\t" << endUm - staUm
         << std::endl;
    ff1 << "unordered_map_time: " << endUm - staUm << " "
        << "mem used: " << endUsedMemUm - startUsedMemUm << "\n";

    double virt = 0.0;
    double res = 0.0;
    for (auto it = configs.begin(); it != configs.end(); it++) {
        // ass is associativity
        size_t ass = it->first;
        size_t bn = it->second;

        staTm = get_usec();

        myTrie::htrie_map<char, uint32_t> hm(ass, bn);
        std::fstream f1(testing_dataset);

        uint64_t startUsedMemTm = getLftMem();
        myTrie::debuging::clear_process_mem_usage();
        while (f1 >> url >> v) {
            hm.insertKV(url, v);
        }
        myTrie::debuging::process_mem_usage(virt, res);

        myTrie::debuging::clear_num();
#ifdef SHRINK_TEST_GROWCUCKOOHASH
        hm.shrink();
        cout << "shrinking cost time: " << shrink_total_time << endl;
#endif
        myTrie::debuging::print_tree_construct<char, uint32_t>(hm.t_root);
        myTrie::debuging::print_tree_construct_v2k<char, uint32_t>(hm.v2k);
        double mem_cal_inside = myTrie::debuging::print_res<char, uint32_t>();

        uint64_t endUsedMemTm = getLftMem();
        f1.close();

        endTm = get_usec();

        cout << "finish trie_map constructing\n";
        cout << "constructing time: " << endTm - staTm << endl;
#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
        cout << "expand cost time: " << expand_cost_time << endl;
#endif

#ifndef TEST_HAT
        cout << "rehash cost time: " << rehash_cost_time << endl;
#endif

        int64_t um_hm_k = 0;
        double percent_k = 0.0;
        double max_percent_k = 0.0;
        double min_percent_k = 0.0;

        // checking:
        for (auto it = m1.begin(); it != m1.end(); it++) {
            std::string url1 = it->first;
            it++;
            std::string url2 = it->first;
            it++;
            std::string url3 = it->first;
            it++;
            std::string url4 = it->first;
            it++;
            std::string url5 = it->first;
            it++;
            if (it == m1.end()) {
                break;
            }

            int64_t hm_get_start = get_usec();
            uint32_t gotfromhm;
            gotfromhm = hm.searchByKey(url1);
            gotfromhm = hm.searchByKey(url2);
            gotfromhm = hm.searchByKey(url3);
            gotfromhm = hm.searchByKey(url4);
            gotfromhm = hm.searchByKey(url5);
            int64_t hm_get_end = get_usec();

            int64_t um_get_start = get_usec();
            uint32_t gotfromum;
            gotfromum = m1[url1];
            gotfromum = m1[url2];
            gotfromum = m1[url3];
            gotfromum = m1[url4];
            gotfromum = m1[url5];
            int64_t um_get_end = get_usec();

            int64_t um_used_time = um_get_end - um_get_start;
            int64_t hm_used_time = hm_get_end - hm_get_start;

            um_hm_k += hm_used_time - um_used_time;

            if (um_used_time == 0) {
                um_used_time = 1;
            }

            double cur_percent_k =
                (double)hm_used_time / (double)um_used_time - 1.0;
            percent_k += cur_percent_k;

            if (max_percent_k < cur_percent_k) {
                max_percent_k = cur_percent_k;
            }
            if (min_percent_k > cur_percent_k) {
                min_percent_k = cur_percent_k;
            }
        }
        cout << "compare to unordered_map: accessing cost diff: "
             << um_hm_k / count << endl;
        ;

        // checking:
        int64_t um_hm_v = 0;
        double percent_v = 0.0;
        double max_percent_v = 0.0;
        double min_percent_v = 0.0;
        for (auto it = m2.begin(); it != m2.end(); it++) {
            uint32_t v1 = it->first;
            it++;
            uint32_t v2 = it->first;
            it++;
            uint32_t v3 = it->first;
            it++;
            uint32_t v4 = it->first;
            it++;
            uint32_t v5 = it->first;
            it++;
            if (it == m2.end()) {
                break;
            }

            int64_t hm_get_start = get_usec();
            std::string gotfromhm = hm.searchByValue(v1);
            gotfromhm = hm.searchByValue(v2);
            gotfromhm = hm.searchByValue(v3);
            gotfromhm = hm.searchByValue(v4);
            gotfromhm = hm.searchByValue(v5);
            int64_t hm_get_end = get_usec();

            int64_t um_get_start = get_usec();
            std::string gotfromum = m2[v1];
            gotfromum = m2[v2];
            gotfromum = m2[v3];
            gotfromum = m2[v4];
            gotfromum = m2[v5];
            int64_t um_get_end = get_usec();

            int64_t um_used_time = um_get_end - um_get_start;
            int64_t hm_used_time = hm_get_end - hm_get_start;

            um_hm_v += hm_used_time - um_used_time;
            if (um_used_time == 0) {
                um_used_time = 1;
            }

            double cur_percent_v =
                (double)hm_used_time / (double)um_used_time - 1.0;
            percent_v += cur_percent_v;
            if (max_percent_v < cur_percent_v) {
                max_percent_v = cur_percent_v;
            }
            if (min_percent_v > cur_percent_v) {
                min_percent_v = cur_percent_v;
            }
        }
        cout << "compare to unordered_map: accessing cost diff: "
             << um_hm_v / count << endl;
        ;

#ifdef TEST_CUCKOOHASH
        ff1 << "cuckoo_hash,,";
#endif
#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
        ff1 << "grow_cuckoo_hash,";
#ifdef REHASH_BEFORE_EXPAND
        ff1 << "rehash_before_expand,";
#else
        ff1 << "expand_before_rehash,";
#endif
#endif
#ifdef TEST_HAT
        ff1 << "hat,,";
#endif
#ifdef TEST_CUCKOOHASH
        ff1 << Associativity << "," << Bucket_num << ","
            << (endTm - staTm) / 1000 / (double)1000 << "," << virt << ","
            << res << "," << mem_cal_inside;
        ff1 << "," << um_hm_k << "," << (double)um_hm_k / (double)count << ","
            << percent_k / (count / 5) * 100.0 << "," << max_percent_k * 100.0
            << "," << min_percent_k * 100.0;
        ff1 << "," << um_hm_v << "," << (double)um_hm_v / (double)count << ","
            << percent_v / (count / 5) * 100.0 << "," << max_percent_v * 100.0
            << "," << min_percent_v * 100.0 << ","
            << (double)rehash_total_num / (double)count << ","
            << ((double)myTrie::debuging::total_pass_trie_node_num /
                (double)count) -
                   1
            << ","
            << (double)myTrie::debuging::hashnode_load /
                   (double)myTrie::debuging::h_n / (double)Max_slot_num *
                   (double)100
            << ","
            << (double)myTrie::debuging::hashnode_max_load /
                   (double)Max_slot_num * (double)100
            << ","
            << (double)myTrie::debuging::hashnode_min_load /
                   (double)Max_slot_num * (double)100
            << "," << 0 << ","
            << (double)rehash_cost_time / (double)1000 / (double)1000
            << ","

#endif

#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
            ff1
            << Associativity << "," << Bucket_num << ","
            << (endTm - staTm) / 1000 / (double)1000 << "," << virt << ","
            << res << "," << mem_cal_inside;
        ff1 << "," << um_hm_k << "," << (double)um_hm_k / (double)count << ","
            << percent_k / (count / 5) * 100.0 << "," << max_percent_k * 100.0
            << "," << min_percent_k * 100.0;
        ff1 << "," << um_hm_v << "," << (double)um_hm_v / (double)count << ","
            << percent_v / (count / 5) * 100.0 << "," << max_percent_v * 100.0
            << "," << min_percent_v * 100.0 << ","
            << (double)rehash_total_num / (double)count << ","
            << ((double)myTrie::debuging::total_pass_trie_node_num /
                (double)count) -
                   1
            << ","
            << (double)myTrie::debuging::hashnode_load /
                   (double)myTrie::debuging::hashnode_total_slot_num *
                   (double)100
            << ","
            << (double)myTrie::debuging::hashnode_max_load /
                   (double)Max_slot_num * (double)100
            << ","
            << (double)myTrie::debuging::hashnode_min_load /
                   (double)Max_slot_num * (double)100
            << "," << (double)expand_cost_time / (double)1000 / (double)1000
            << "," << (double)rehash_cost_time / (double)1000 / (double)1000
            << ","

#endif
#ifdef TEST_HAT
            ff1
            << hm.burst_threshold << "," << hm.bucket_num << ","
            << (endTm - staTm) / 1000 / (double)1000 << "," << virt << ","
            << res << "," << mem_cal_inside;
        ff1 << "," << um_hm_k << "," << um_hm_k / count << ","
            << percent_k / (count / 5) * 100.0 << "," << max_percent_k * 100.0
            << "," << min_percent_k * 100.0;
        ff1 << "," << um_hm_v << "," << um_hm_v / count << ","
            << percent_v / (count / 5) * 100.0 << "," << max_percent_v * 100.0
            << "," << min_percent_v * 100.0
            << ","
            // << (double)rehash_total_num / (double)count << ","
            << ((double)myTrie::debuging::total_pass_trie_node_num /
                (double)count) -
                   1
            << ","
            << (double)myTrie::debuging::hashnode_load /
                   (double)myTrie::debuging::h_n / (double)hm.burst_threshold
            << ","
#endif
            << myTrie::debuging::t_n << "," << myTrie::debuging::h_n << endl;

        ff1.flush();

#ifdef TEST_CUCKOOHASH
        rehash_cost_time = 0;
        rehash_total_num = 0;
#endif

#ifdef IMPROVE_BURST
        fstream recal_time_file("recal_record", ios::out | ios::app);
        recal_time_file << testing_dataset << "," << Associativity << ","
                        << Bucket_num << ",";
        recal_time_file << recal_element_num_of_1st_char_counter;
        recal_time_file << "," << burst_total_counter << endl;
        recal_time_file.flush();
        recal_element_num_of_1st_char_counter = 0;
        burst_total_counter = 0;
#endif

#if (defined SHRINK_TEST_GROWCUCKOOHASH) || (defined TEST_GROWCUCKOOHASH)
        expand_cost_time = 0;
        rehash_cost_time = 0;
        rehash_total_num = 0;
#endif
#ifdef SHRINK_TEST_GROWCUCKOOHASH
        shrink_total_time = 0;
#endif

        // correctness check
        if (test_and_print_wrong_test) {
            vector<string> wrong_search_key;
            vector<uint32_t> wrong_search_value;

            for (auto it = m1.begin(); it != m1.end(); it++) {
                if (it->second != hm.searchByKey(it->first)) {
                    wrong_search_key.push_back(it->first);
                }
            }

            cout << "test key finish!\n";
            cout << "wrong key_searching num: " << wrong_search_key.size()
                 << endl;

            for (auto it = m2.begin(); it != m2.end(); it++) {
                if (it->second != hm.searchByValue(it->first)) {
                    wrong_search_value.push_back(it->first);
                }
            }

            cout << "test value finish!\n";
            cout << "wrong value_searching num: " << wrong_search_value.size()
                 << endl;

            if (wrong_search_value.size() != 0 ||
                wrong_search_key.size() != 0) {
                cout << "testing failed\n";
                exit(0);
            }
        }

        if (manually_test) {
            while (cin >> url) {
                cout << "get value: " << hm.searchByKey(url) << endl;
            }
        }
    }
    ff1.close();
}
