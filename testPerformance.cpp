#include <boost/unordered_map.hpp>
#include <map>
#include <vector>
// #include "myTrie.hpp"
// debug
#include <fstream>
#include "debug.hpp"
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
    elem_per_buck.push_back(11);
    elem_per_buck.push_back(31);
    elem_per_buck.push_back(101);

    // analyse by bucket_num
    bucket_nums.push_back(5);
    bucket_nums.push_back(31);
    bucket_nums.push_back(53);

    bucket_nums.push_back(103);
    bucket_nums.push_back(991);
    bucket_nums.push_back(1999);
    bucket_nums.push_back(4139);

    // configs: pair<burst_point, bucket_num>
    for (auto it = elem_per_buck.begin(); it != elem_per_buck.end(); it++) {
        for (auto itt = bucket_nums.begin(); itt != bucket_nums.end(); itt++) {
            configs.push_back(std::pair<size_t, size_t>((*it * *itt), *itt));
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

    std::fstream f("dataset/str_normal");
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
        size_t bp = it->first;
        size_t bn = it->second;

        staTm = get_usec();

        myTrie::htrie_map<char, uint32_t> hm(bp, bn);
        std::fstream f1("dataset/str_normal");

        uint64_t startUsedMemTm = getLftMem();
        myTrie::debuging::clear_process_mem_usage();
        while (f1 >> url >> v) {
            hm.insertKV(url, v);
        }
        myTrie::debuging::process_mem_usage(virt, res);

        myTrie::debuging::print_tree_construct<char, uint32_t>(hm.t_root);
        double mem_cal_inside = myTrie::debuging::print_res<char, uint32_t>();

        uint64_t endUsedMemTm = getLftMem();
        f1.close();

        endTm = get_usec();

        cout << "finish trie_map constructing\n";

        int64_t um_hm_k = 0;
        double percent_k = 0.0;
        double max_percent_k = 0.0;
        double min_percent_k = 0.0;

        // checking:
        std::fstream f3("test_res_wrong_key", std::ios::out | std::ios::app);
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

            // correctness check
            // if (gotfromhm == it->second) {
            // } else {
            //     f3 << "key_wrong answer: " << it->first << " got "
            //        << hm.searchByKey(it->first).second
            //        << " from hm , actual value: " << it->second << std::endl;
            //     uint32_t v = it->second;
            //     f3 << "got from hm: " << hm.searchByKey(it->first).second;
            //     f3 << "\ngot from file: " << v;
            //     for (size_t i = 0; i != sizeof(v); i++) {
            //         f3 << (unsigned int)*((char*)(&v + sizeof(char) * i))
            //            << ",";
            //     }
            //     f3 << std::endl;
            //     f3.flush();
            //     std::cout << "wrong answer!\n";
            //     exit(0);
            // }
        }
        f3.close();
        cout << "compare to unordered_map: accessing cost diff: "
             << um_hm_k / count << endl;
        ;

        // checking:
        int64_t um_hm_v = 0;
        double percent_v = 0.0;
        double max_percent_v = 0.0;
        double min_percent_v = 0.0;
        std::fstream f4("test_res_wrong_value", std::ios::out | std::ios::app);
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

            double cur_percent_v = hm_used_time / um_used_time - 1;
            percent_v += cur_percent_v;
            if (max_percent_v < cur_percent_v) {
                max_percent_v = cur_percent_v;
            }
            if (min_percent_v > cur_percent_v) {
                min_percent_v = cur_percent_v;
            }

            // correctness check
            // if (gotfromhm == it->second) {
            // } else {
            //     f4 << "value_wrong answer: " << it->first << " got "
            //        << hm.searchByValue(it->first)
            //        << " from hm , actual value: " << it->second << std::endl;
            //     std::string v = it->second;
            //     f4 << "got from hm: " << hm.searchByValue(it->first);
            //     f4 << "\ngot from file: " << v;
            //     f4 << std::endl;
            //     f4.flush();
            //     std::cout << "wrong answer!\n";
            //     exit(0);
            // }
        }
        f4.close();
        cout << "compare to unordered_map: accessing cost diff: "
             << um_hm_v / count << endl;
        ;

        // ff1 << hm.burst_threshold << "," << hm.bucket_num << "\t init: ,"
        //     << endTm - staTm << ","
        //     << "virt : ," << virt << ",\t"
        //     << "res : ," << res << ",\t"
        //     << "inside_mem: ," << mem_cal_inside << ",\t";
        // ff1 << "key total: ," << um_hm_k << ",\t avg: ," << um_hm_k / count
        //     << ",\t perc: ," << percent_k / (count / 5) * 100.0 << ","
        //     << " max_%: ," << max_percent_k << ", min_%: ," << min_percent_k;
        // ff1 << " value total: ," << um_hm_v << ",\t avg: ," << um_hm_v /
        // count
        //     << ",\t perc: ," << percent_v / (count / 5) * 100.0 << ","
        //     << " max_%: ," << max_percent_v << ", min_%: ," << min_percent_v
        //     << endl
        //     << endl;

        ff1 << bp / bn << "," << bn << "," << endTm - staTm << "," << virt
            << "," << res << "," << mem_cal_inside;
        ff1 << "," << um_hm_k << "," << um_hm_k / count << ","
            << percent_k / (count / 5) * 100.0 << "," << max_percent_k * 100.0
            << "," << min_percent_k * 100.0;
        ff1 << "," << um_hm_v << "," << um_hm_v / count << ","
            << percent_v / (count / 5) * 100.0 << "," << max_percent_v * 100.0
            << "," << min_percent_v * 100.0 << endl;

        //--------------printing wrong result-------------------
        std::ifstream f4wrong1("test_res_wrong_key", std::ios::in);
        char line[1024] = {0};
        while (f4wrong1.getline(line, sizeof(line))) {
            std::cerr << line << std::endl;
        }

        std::ifstream f4wrong2("test_res_wrong_value", std::ios::in);
        while (f4wrong2.getline(line, sizeof(line))) {
            std::cerr << line << std::endl;
        }

        startUsedMemTm = getLftMem();
        // hm.deleteMyself();
        endUsedMemTm = getLftMem();
        cout << "finish checking and printed correct/wrong "
                "res\n--------------------------------------\n";
        // }
        // while(cin >> url){
        //     cout << "get value: " << hm.searchByKey(url) << endl;
        // }
    }

    ff1.close();
}
