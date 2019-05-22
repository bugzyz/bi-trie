#include "../unified_impl/test_switcher.hpp"
#include <iostream>
using namespace std;

int main() {
    zyz_trie::bi_trie<char, int, 3, 2> bt;

    int cin_time = 0;
    cout << "Input element number\n";
    cin >> cin_time;

    for(int i = 0; i!=cin_time;i++){
        string key = "";
        int value = 0;
        cin >> key >> value;
        bt.insert_kv(key, value);
    }

    cout << "Input finished\n";

    for(int i=0;i!=cin_time; i++) {
        string key = "";
        cin >> key;
        cout << bt[key]<< endl;
    }

    cout << "Test search by key finished\n";


    for(int i=0;i!=cin_time; i++) {
        int value = 0;
        cin >> value;
        cout << bt[value]<< endl;
    }

    cout << "Test search by value finished\n";
}