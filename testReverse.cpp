#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "myTrie.hpp"

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
    myTrie::htrie_map<char, uint32_t> hm(16384, 8192);
    map<string, uint32_t> m1;
    map<uint32_t, string> m2;
    fstream f("dataset/str_normal");
    string url;
    uint32_t v;
    uint64_t sta = get_usec();
    while (f >> url >> v) {
        hm.insertKV(url, v);
        m1[url] = v;
        m2[v] = url;
    }
    uint64_t end = get_usec();
    std::cout << "loaded 3 map in " << (end - sta) / 1000 << " ms" << std::endl;

    uint64_t count = m1.size();
    uint64_t searchByK_count_good = 0;
    uint64_t searchByK_count_bad = 0;
    uint64_t searchByV_count_good = 0;
    uint64_t searchByV_count_bad = 0;

    for (auto it = m1.begin(); it != m1.end(); it++) {
        if (hm.searchByKey(it->first) == it->second) {
            std::cout << "good\n";
            std::cout << "ans: " << it->second << std::endl;
            std::cout << "got " << hm.searchByKey(it->first) << std::endl;
            searchByK_count_good++;
        } else {
            std::cout << "wrong\n";
            std::cout << "ans: " << it->second << std::endl;
            std::cout << "got " << hm.searchByKey(it->first) << std::endl;
            searchByK_count_bad++;
        }
    }
    std::cout << "-------------------\n";

    for (auto it = m2.begin(); it != m2.end(); it++) {
        if (hm.searchByValue(it->first) == it->second) {
            std::cout << "good\n";
            std::cout << "ans: " << it->second << std::endl;
            std::cout << "got " << hm.searchByValue(it->first) << std::endl;
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