bool test_and_print_wrong_test = true;
bool manually_test = false;

// test which
#include "unified_impl/5_prototype_cuckoo_grow_buc_shrink.hpp"
// #include "unified_impl/4_prototype_cuckoo_grow_ass_shrink.hpp"
// #include "unified_impl/3_prototype_cuckoo_shrink.hpp"
// #include "unified_impl/2_prototype_shrink.hpp"
// #include "unified_impl/1_tessil_hat_impl.hpp"

// debug
#include <fstream>

#include <boost/unordered_map.hpp>
#include <map>
#include <set>
#include <vector>

#include <stdint.h>

#include <fstream>
#include <iostream>

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

    // configs: pair<bucket_num, elem_per_bucket>
    // stable bucket_num
    for (auto itt = elem_per_buck.begin(); itt != elem_per_buck.end(); itt++) {
        for (auto it = bucket_nums.begin(); it != bucket_nums.end(); it++) {
            configs.push_back(std::pair<size_t, size_t>(*itt, *it));
        }
    }

    set<string> datasets;
    datasets.insert("dataset/str_normal");
    datasets.insert("dataset/id_yago/cut_str_normal");

    for (auto set_it = datasets.begin(); set_it != datasets.end(); set_it++) {
        string testing_dataset = *set_it;

        cout << "testing file: " << testing_dataset << endl;
        //------------testing--------------
        std::string url;
        uint32_t v;

        std::fstream f(testing_dataset);
        boost::unordered_map<string, uint32_t> m1;
        boost::unordered_map<uint32_t, string> m2;

        while (f >> url >> v) {
            m1[url] = v;
            m2[v] = url;
        }
        f.close();

        for (auto it = configs.begin(); it != configs.end(); it++) {
            // ass is associativity
            size_t ass = it->first;
            size_t bn = it->second;

            cout << "running config(" << ass << "," << bn << ")\n";

            myTrie::htrie_map<char, uint32_t> hm(ass, bn);

            uint64_t start = get_time();
            for (auto itt = m1.begin(); itt != m1.end(); itt++) {
                hm.insertKV(itt->first, itt->second);
            }
            uint64_t end = get_time();

            cout << "trie-map insertion finished! tooks "
                 << ((double)((end - start) / 1000)) / (double)1000 << " s\n";

            hm.shrink();

            // correctness check
            if (test_and_print_wrong_test) {
                vector<string> wrong_search_key;
                vector<uint32_t> wrong_search_value;

                size_t counter1 = 0;
                for (auto test_it = m1.begin(); test_it != m1.end();
                     test_it++) {
                    counter1++;
                    if (test_it->second != hm.searchByKey(test_it->first)) {
                        wrong_search_key.push_back(test_it->first);
                    }
                }

                cout << "test key finish!\n";
                cout << "wrong key_searching num: " << wrong_search_key.size()
                     << "/" << counter1 << endl;

                size_t counter2 = 0;
                for (auto test_it = m2.begin(); test_it != m2.end();
                     test_it++) {
                    counter2++;
                    if (test_it->second != hm.searchByValue(test_it->first)) {
                        wrong_search_value.push_back(test_it->first);
                    }
                }

                cout << "test value finish!\n";
                cout << "wrong val_searching num: " << wrong_search_value.size()
                     << "/" << counter2 << endl;

                if (wrong_search_value.size() != 0 ||
                    wrong_search_key.size() != 0) {
                    cout << "testing failed\n";
                    exit(0);
                }
            }

            cout << "-----------------\n";
        }
    }
}
