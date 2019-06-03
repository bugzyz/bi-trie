# Bi-trie
The bi-trie stands for bidirectional mapping trie.
Bi-trie provides bidirectional mapping between key(string) and value(int, uint32_t, unsigned int ...) in a very memory-efficient way.

This implementation is based on a burst trie design that a leaf node is a hash table to store elements that greatly reduces the number of intermediate node and pointer traversal in a trie searching traversal.
