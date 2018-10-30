#include <iostream>
#include <vector>
#include "myTrie.hpp"

namespace myTrie {
namespace debuging {
using std::cout;
using std::cout;
using std::endl;
using std::endl;
uint64_t h_n;
uint64_t t_n;
template <typename CharT, typename T>
void print_tree_construct(typename myTrie::htrie_map<CharT, T>::anode* root) {
    if (root == nullptr) {
        return;
    }
    // print bucket
    if (root->isHashNode()) {
        h_n++;
    } else if (root->isTrieNode()) {
        std::map<CharT, typename myTrie::htrie_map<CharT, T>::anode*> childs =
            ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)->childs;
        if (childs.size() == 0) {
            print_tree_construct(
                ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                    ->getOnlyHashNode());
        } else {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                t_n++;
                print_tree_construct(it->second);
            }
        }
    } else {
        std::cout << "node is not trie nor hash node\n";
        exit(0);
    }
    return;
}

void print_res() {
    std::cout << "trie_node: " << t_n << " hash_node: " << h_n << std::endl;
    t_n = 0;
    h_n = 0;
}
}  // namespace debuging
}  // namespace myTrie