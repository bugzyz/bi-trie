#include "myTrie.hpp"
#include <vector>
#include <iostream>

namespace myTrie {  //------------------------DEBUG, CORRECTNESS
                    //CHECK------------------------------
namespace debuging {
using std::cout;
using std::endl;
template <class CharT, class T>
void print_tree_struct(typename myTrie::htrie_map<CharT, T>::anode* root) {
    debugStream << "in node: " << (void*)root
                << " with type: " << (root->isHashNode() ? "hash" : "trie")
                << std::endl;
    // print bucket
    if (root->isHashNode()) {
        debugStream << "searching in hashnode" << endl;
        std::map<size_t, std::vector<std::pair<std::string, T>>> buckets_info =
            ((typename myTrie::htrie_map<CharT, T>::hash_node*)root)
                ->get_hash_nodes_buckets_info();
        for (auto it = buckets_info.begin(); it != buckets_info.end(); it++) {
            debugStream << it->first << ":";
            std::vector<std::pair<std::string, T>> bucket_string = it->second;
            for (int i = 0; i != bucket_string.size(); i++) {
                debugStream << bucket_string[i].first << ",";
            }
            debugStream << "|" << endl;
        }
    } else if (root->isTrieNode()) {
        std::map<CharT, typename myTrie::htrie_map<CharT, T>::anode*> childs =
            ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                ->getChildsMap();
        for (auto it = childs.begin(); it != childs.end(); it++) {
            debugStream << it->first << "\t\t";
        }
        debugStream << "\nnext level:\n";
        if (childs.size() == 0) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_struct<CharT, T>(it->second);
            }
        } else {
            debugStream << "size of child is 0" << endl;
            debugStream << "fake_root is " << (void*)root << endl;
            // print_tree_struct<CharT, T>(
            //     ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
            //         ->getOnlyHashNode());
        }
    } else {
        debugStream << "node is not trie nor hash node\n";
        exit(0);
    }
}

uint32_t trie_node_num;
uint32_t hash_node_num;
uint64_t bucket_buf_size;

template <class CharT, class T>
void print_tree_use_mem(typename myTrie::htrie_map<CharT, T>::anode* root) {
    // print bucket
    if (root->isHashNode()) {
        debugStream << "searching in hashnode" << endl;
        std::map<size_t, std::vector<std::pair<std::string, T>>> buckets_info =
            ((typename myTrie::htrie_map<CharT, T>::hash_node*)root)
                ->get_hash_nodes_buckets_info();
        for (auto it = buckets_info.begin(); it != buckets_info.end(); it++) {
            debugStream << it->first << ":";
            std::vector<std::pair<std::string, T>> bucket_string = it->second;
            for (int i = 0; i != bucket_string.size(); i++) {
                debugStream << bucket_string[i].first << ",";
            }
            debugStream << "|" << endl;
        }
    } else if (root->isTrieNode()) {
        std::map<CharT, typename myTrie::htrie_map<CharT, T>::anode*> childs =
            ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
                ->getChildsMap();
        for (auto it = childs.begin(); it != childs.end(); it++) {
            debugStream << it->first << "\t\t";
        }
        debugStream << "\nnext level:\n";
        if (childs.size() == 0) {
            for (auto it = childs.begin(); it != childs.end(); it++) {
                print_tree_struct<CharT, T>(it->second);
            }
        } else {
            debugStream << "size of child is 0" << endl;
            debugStream << "fake_root is " << (void*)root << endl;
            // print_tree_struct<CharT, T>(
            //     ((typename myTrie::htrie_map<CharT, T>::trie_node*)root)
            //         ->getOnlyHashNode());
        }
    } else {
        debugStream << "node is not trie nor hash node\n";
        exit(0);
    }
}

template <class CharT, class T>
void print_htrie_map(htrie_map<CharT, T> hm,
                     std::vector<std::string> checklist) {
    // correctness check
    debugStream << "finding func check:\n";
    for (auto it = checklist.begin(); it != checklist.end(); it++) {
        debugStream << "search " << *it << " got " << hm[*it] << endl;
    }

    debugStream << "-------------------------\n";

    // node check
    // print_tree_struct<CharT, T>(hm.t_root);
}

template <class CharT>
void printDiff(const CharT* key_lhs, std::size_t key_size_lhs,
               const CharT* key_rhs, std::size_t key_size_rhs, bool equal) {
    debugStream << "==comparing: \n";
    for (size_t i = 0; i != key_size_lhs; i++) {
        // debugStream << (unsigned int)*(key_lhs + i) << ',';
        debugStream << *(key_lhs + i) << ',';
    }
    debugStream << "with size:" << key_size_lhs << "\n and \n";
    for (size_t i = 0; i != key_size_rhs; i++) {
        // debugStream << (unsigned int)*(key_rhs + i) << ',';
        debugStream << *(key_rhs + i) << ',';
    }
    debugStream << "with size:" << key_size_rhs << std::endl;
    debugStream << " keyEqual res: "
                << (std::memcmp(key_lhs, key_rhs,
                                key_size_lhs * sizeof(CharT)) == 0)
                << " keyEuqal res: " << equal << std::endl;
}
}  // namespace debuging
}  // namespace myTrie