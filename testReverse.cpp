#include <fstream>
#include <iostream>
#include <string>
#include "myTrie.hpp"
using namespace std;
int main() {
  myTrie::htrie_map<char, uint32_t> hm(80,5);
  fstream f("str_normal.txt");
  string url;
  uint32_t v;
  while (f >> url >> v) {
      hm.insertKV(url,v);
  }
}