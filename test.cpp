#include <sys/sysinfo.h>
#include <boost/unordered_map.hpp>
#include <map>
#include <vector>
#include "myTrie.hpp"
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

uint64_t getLftMem() {
    struct sysinfo inf;
    sysinfo(&inf);

    return inf.mem_unit *
           (inf.totalram + inf.totalswap - inf.freeram - inf.freeswap);
    // return get_proc_mem(getpid())+get_proc_virtualmem(getpid());
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
    std::cout << "starto!\n";
    // test data
    vector<uint32_t> bursts;
    vector<uint32_t> ratios;

    // bursts.push_back(80);
    // bursts.push_back(150);
    // bursts.push_back(300);
    // bursts.push_back(550);
    // bursts.push_back(800);
    // bursts.push_back(999);
    // bursts.push_back(9800);

    //todo
    // bursts.push_back(1500);
    // bursts.push_back(4500);
    // bursts.push_back(7200);
    // bursts.push_back(11001);
    // bursts.push_back(13200);
    //
    // bursts.push_back(16384);
    // bursts.push_back(20000);
    // bursts.push_back(25000);
    bursts.push_back(36000);
    bursts.push_back(42000);
    bursts.push_back(50000);

    ratios.push_back(1);
    ratios.push_back(2);
    ratios.push_back(3);
    ratios.push_back(4);
    ratios.push_back(6);
    ratios.push_back(8);
    ratios.push_back(11);
    ratios.push_back(21);
    ratios.push_back(37);
    ratios.push_back(50);
    ratios.push_back(100);
    ratios.push_back(300);
    ratios.push_back(501);
    ratios.push_back(750);
    ratios.push_back(1000);
    ratios.push_back(2000);
    ratios.push_back(0);

    fstream ff1("result_key", std::ios::out);
    fstream ff2("result_value", std::ios::out);
    // burst_threshold
    vector<unsigned int> bt;
    // bucket_num
    vector<unsigned int> bn;

    bt.push_back(80);

    //------------testing--------------
    std::string url;
    uint32_t v;
    uint64_t staUm;
    uint64_t endUm;
    uint64_t staTm;
    uint64_t endTm;

    std::fstream f("str_normal");
    boost::unordered_map<string, uint32_t> m;
    boost::unordered_map<uint32_t, string> m2;
    // std::map<string, uint32_t> m;

    staUm = get_usec();
    uint64_t startUsedMemUm = getLftMem();

    uint32_t count = 0;
    while (f >> url >> v) {
        m[url] = v;
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
    ff2 << "unordered_map_time: " << endUm - staUm << " "
        << "mem used: " << endUsedMemUm - startUsedMemUm << "\n";
    uint32_t last_time_bn = 0;
    for (int i = 0; i != bursts.size(); i++) {
        for (int j = 0; j != ratios.size(); j++) {
            uint32_t bp = bursts[i];
            uint32_t bn = 0;
            if (ratios[j] == 0)
                bn = 1;
            else
                bn = bp / ratios[j];
            if (bn == 0) continue;

            if (bn == last_time_bn) {
                continue;
            }

            last_time_bn = bn;

            staTm = get_usec();
            // todo:
            myTrie::htrie_map<char, uint32_t> hm(bp, bn);
            std::fstream f1("str_normal");

            uint64_t startUsedMemTm = getLftMem();
            while (f1 >> url >> v) {
                hm.insertKV(url, v);
            }
            uint64_t endUsedMemTm = getLftMem();
            f1.close();

            endTm = get_usec();

            cout << "myTrie(" << hm.burst_threshold << ",\t " << hm.bucket_num
                 << ") use time: usec: \t\t" << endTm - staTm << std::endl;
            cout << "finish trie_map constructing\n";

            int64_t um_hm = 0;
            // checking:
            std::fstream f3("test_res_wrong_key",
                            std::ios::out | std::ios::app);
            for (auto it = m.begin(); it != m.end(); it++) {
                int64_t hm_get_start = get_usec();
                uint32_t gotfromhm = hm.searchByKey(it->first);
                int64_t hm_get_end = get_usec();

                int64_t um_get_start = get_usec();
                uint32_t gotfromum = m[it->first];
                int64_t um_get_end = get_usec();

                um_hm +=
                    (um_get_end - um_get_start) - (hm_get_end - hm_get_start);

                if (gotfromhm == it->second) {
                } else {
                    f3 << "key_wrong answer: " << it->first << " got "
                       << hm[it->first]
                       << " from hm , actual value: " << it->second
                       << std::endl;
                    uint32_t v = it->second;
                    f3 << "got from hm: " << hm.searchByKey(it->first);
                    f3 << "\ngot from file: " << v;
                    for (size_t i = 0; i != sizeof(v); i++) {
                        f3 << (unsigned int)*((char*)(&v + sizeof(char) * i))
                           << ",";
                    }
                    f3 << std::endl;
                    f3.flush();
                    std::cout << "wrong answer!\n";
                    exit(0);
                }
            }
            f3.close();
            cout << "compare to unordered_map: accessing cost diff: "
                 << um_hm / count << endl;
            ;
            ff1 << "search by key: " << hm.burst_threshold << ","
                << hm.bucket_num << "\t use time: usec: ," << endTm - staTm
                << ","
                << "mem used: ," << endUsedMemTm - startUsedMemTm << ",\t";
            ff1 << "accessing cost total diff: ," << um_hm
                << ",\t avg: " << um_hm / count << endl;
            ;

            // checking:
            std::fstream f4("test_res_wrong_value",
                            std::ios::out | std::ios::app);
            for (auto it = m2.begin(); it != m2.end(); it++) {
                int64_t hm_get_start = get_usec();
                std::string gotfromhm = hm.searchByValue(it->first);
                int64_t hm_get_end = get_usec();

                int64_t um_get_start = get_usec();
                std::string gotfromum = m2[it->first];
                int64_t um_get_end = get_usec();

                um_hm +=
                    (um_get_end - um_get_start) - (hm_get_end - hm_get_start);

                if (gotfromhm == it->second) {
                } else {
                    f4 << "value_wrong answer: " << it->first << " got "
                       << hm.searchByValue(it->first)
                       << " from hm , actual value: " << it->second
                       << std::endl;
                    std::string v = it->second;
                    f4 << "got from hm: " << hm.searchByValue(it->first);
                    f4 << "\ngot from file: " << v;
                    f4 << std::endl;
                    f4.flush();
                    std::cout << "wrong answer!\n";
                    exit(0);
                }
            }
            f4.close();
            cout << "compare to unordered_map: accessing cost diff: "
                 << um_hm / count << endl;
            ;
            ff2 << "search by value: " << hm.burst_threshold << ","
                << hm.bucket_num << "\t use time: usec: ," << endTm - staTm
                << ","
                << "mem used: ," << endUsedMemTm - startUsedMemTm << ",\t";
            ff2 << "accessing cost total diff: ," << um_hm
                << ",\t avg: " << um_hm / count << endl;
            ;
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

            hm.deleteMyself();
            endUsedMemTm = getLftMem();
            cout << "release: " << startUsedMemTm - endUsedMemTm << "\n";
            cout << "finish checking and printed correct/wrong "
                    "res\n--------------------------------------\n";
        }
    }

    ff1.close();
    ff2.close();

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
