#include "myTrie.hpp"
#include <boost/unordered_map.hpp>
#include <map>
#include <vector>

// debug
#include <fstream>
using std::fstream;

#include <stdint.h>
#include <sys/time.h>
#include <unistd.h>
static uint64_t get_usec() {
    struct timespec tp;
    /* POSIX.1-2008: Applications should use the clock_gettime() function
       instead of the obsolescent gettimeofday() function. */
    /* NOTE: The clock_gettime() function is only available on Linux.
       The mach_absolute_time() function is an alternative on OSX. */
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return ((tp.tv_sec * 1000 * 1000) + (tp.tv_nsec / 1000));
}

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
    fstream ff("result", std::ios::out);
    // burst_threshold
    vector<unsigned int> bt;
    // bucket_num
    vector<unsigned int> bn;

    bt.push_back(80);

    //------------testing--------------
    std::string url;
    uint32_t v;
    uint64_t sta;
    uint64_t end;

    std::fstream f("str_normal");
    boost::unordered_map<string, uint32_t> m;
    boost::unordered_map<uint32_t, string> m2;
    // std::map<string, uint32_t> m;

    sta = get_usec();
    uint32_t count = 0;
    while (f >> url >> v) {
        m[url] = v;
        // todo
        m2[v] = url;
        count++;
    }
    end = get_usec();
    f.close();
    cout << "unordered_map use time: usec: \t\t\t" << end - sta << std::endl;
    ff << "unordered_map use time: usec: \t\t\t" << end - sta << std::endl;

    for (int i = 20000; i != 100000; i += 20000) {
        sta = get_usec();
        // todo:
        myTrie::htrie_map<char, uint32_t> hm(i, i / 400);
        std::fstream f1("str_normal");

        while (f1 >> url >> v) {
            // cout << "---------------------------------------> working on
            // num."
            //      << count << " url: " << url << " set value: " << v << endl;
            // hm[url] = v;
            hm.insertKV(url, v);
        }
        f1.close();

        end = get_usec();
        ff << "myTrie(" << hm.burst_threshold << ",\t " << hm.bucket_num
           << ") use time: usec: \t\t" << end - sta << std::endl;
        cout << "myTrie(" << hm.burst_threshold << ",\t " << hm.bucket_num
             << ") use time: usec: \t\t" << end - sta << std::endl;
        cout << "finish trie_map constructing\n";

        m.begin()->second = 123456;

        int64_t um_hm = 0;
        // checking:
        std::fstream f3("test_res_wrong", std::ios::out);
        for (auto it = m.begin(); it != m.end(); it++) {
            int64_t hm_get_start = get_usec();
            uint32_t gotfromhm = hm.searchByKey[it->first];
            int64_t hm_get_end = get_usec();

            int64_t um_get_start = get_usec();
            uint32_t gotfromum = m[it->first];
            int64_t um_get_end = get_usec();

            um_hm += (um_get_end - um_get_start) - (hm_get_end - hm_get_start);

            if (gotfromhm == it->second) {
            } else {
                f3 << "wrong answer: " << it->first << " got " << hm[it->first]
                   << " from hm , actual value: " << it->second << std::endl;
                uint32_t v = it->second;
                f3 << "got from hm: ";
                for (size_t i = 0; i != sizeof(v); i++) {
                    f3 << (unsigned int)*(
                              (char*)(&(hm[it->first]) + sizeof(char) * i))
                       << ",";
                }
                f3 << "\ngot from file: ";
                for (size_t i = 0; i != sizeof(v); i++) {
                    f3 << (unsigned int)*((char*)(&v + sizeof(char) * i))
                       << ",";
                }
                f3 << std::endl;
                f3.flush();
            }
        }
        f3.close();
        cout << "compare to unordered_map: accessing cost diff: "
             << um_hm / count << endl;;
        ff << std::endl;

// ------------------------------------------------------------

        um_hm = 0;
        // checking:
        std::fstream f3("test_res_wrong", std::ios::out);
        for (auto it = m2.begin(); it != m2.end(); it++) {
            int64_t hm_get_start = get_usec();
            std::string gotfromhm = hm.searchByValue[it->first];
            int64_t hm_get_end = get_usec();

            int64_t um_get_start = get_usec();
            uint32_t gotfromum = m2[it->first];
            int64_t um_get_end = get_usec();

            um_hm += (um_get_end - um_get_start) - (hm_get_end - hm_get_start);

            if (gotfromhm == it->second) {
            } else {
                f3 << "wrong answer: " << it->first << " got " << hm[it->first]
                   << " from hm , actual value: " << it->second << std::endl;
                uint32_t v = it->second;
                f3 << "got from hm: ";
                for (size_t i = 0; i != sizeof(v); i++) {
                    f3 << (unsigned int)*(
                              (char*)(&(hm[it->first]) + sizeof(char) * i))
                       << ",";
                }
                f3 << "\ngot from file: ";
                for (size_t i = 0; i != sizeof(v); i++) {
                    f3 << (unsigned int)*((char*)(&v + sizeof(char) * i))
                       << ",";
                }
                f3 << std::endl;
                f3.flush();
            }
        }
        f3.close();
        cout << "compare to unordered_map: accessing cost diff: "
             << um_hm / count << endl;;
        ff << std::endl;


        //--------------printing wrong result-------------------
        std::ifstream fofwrong("test_res_wrong", std::ios::in);
        char line[1024] = {0};
        while (fofwrong.getline(line, sizeof(line))) {
            std::cerr << line << std::endl;
        }
        cout << "finish checking and printed correct/wrong "
                "res\n--------------------------------------\n";
    }

    // while (true) {
    //     cout << "check url: ";
    //     cout.rdbuf(fileBufcc);
    //     cin >> url;
    //     sta = get_usec();
    //     cout << "check url: " << url << " got " << hm[url] << endl;
    //     end = get_usec();
    //     cout << "use " << (end - sta) / 1000 << "ms" << endl;
    //     cout.rdbuf(coutBuf);
    // }
}