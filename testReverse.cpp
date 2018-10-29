#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "myTrie.hpp"
using namespace std;
int main() {
    streambuf* coutBuf = cout.rdbuf();
    myTrie::htrie_map<char, uint32_t> hm(50, 5);
    map<string, uint32_t> m1;
    map<uint32_t, string> m2;
    fstream f("str_normal");
    string url;
    uint32_t v;
    while (f >> url >> v) {
        ofstream of("out.txt", std::ios::app);
        // cout << "working on " << url << endl;
        if (v == 132873) {
            streambuf* buf = of.rdbuf();
            cout.rdbuf(buf);
        }
        hm.insertKV(url, v);
        if (v == 132873) {
            cout.rdbuf(coutBuf);
            of.close();
        }
        m1[url] = v;
        m2[v] = url;
    }

    uint64_t count = m1.size();
    uint64_t searchByK_count_good = 0;
    uint64_t searchByK_count_bad = 0;
    uint64_t searchByV_count_good = 0;
    uint64_t searchByV_count_bad = 0;

    for (auto it = m1.begin(); it != m1.end(); it++) {
        if (hm.searchByKey(it->first) == it->second) {
            // std::cout << "good\n";
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
        ofstream of("out.txt", std::ios::app);

        // if (it->first == 132873) {
        //     streambuf* buf = of.rdbuf();

        //     cout.rdbuf(buf);
        // }
        if (hm.searchByValue(it->first) == it->second) {
            // std::cout << "good\n";
            searchByV_count_good++;

        } else {
            std::cout << "wrong\n";
            std::cout << "ans: " << it->second << std::endl;
            std::cout << "got " << hm.searchByValue(it->first) << std::endl;
            searchByV_count_bad++;
        }
        // if (it->first == 132873) {
        //     cout.rdbuf(coutBuf);
        //     of.close();
        // }
    }
    std::cout << "-------------------\n";
    for (auto it = m2.begin(); it != m2.end(); it++) {
        if (hm.searchByValue(it->first) == it->second) {
            // std::cout << "good\n";
            searchByV_count_good++;

        } else {
            std::cout << "wrong\n";
            std::cout << "ans: " << it->second << std::endl;
            std::cout << "got " << hm.searchByValue(it->first) << std::endl;
            searchByV_count_bad++;
        }
    }
    cout << "total: " << count << endl;
    cout << "check key: good:" << searchByK_count_good
         << " bad:" << searchByK_count_bad << std::endl;
    cout << "check key: good:" << searchByV_count_good
         << " bad:" << searchByV_count_bad << std::endl;
    cout << "finish\n";
    // cout << ">";
    // while (cin >> v) {
    //     cout << ">";
    //     hm.searchByValue(v);
    // }
}