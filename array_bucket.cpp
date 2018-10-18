#include <cstdint>
// for std::strlen
#include <cstring>
#include <map>
#include <vector>

#include <string>

// numeric_limit
#include <limits>

#include <iostream>
using std::cout;
using std::endl;

#define BUCKET_INIT_COUNT 5

namespace myTrie {
namespace hashRelative {

template <class CharT>
std::size_t hash(const CharT* key, std::size_t key_size) {
    static const std::size_t init = std::size_t(
        (sizeof(std::size_t) == 8) ? 0xcbf29ce484222325 : 0x811c9dc5);
    static const std::size_t multiplier =
        std::size_t((sizeof(std::size_t) == 8) ? 0x100000001b3 : 0x1000193);

    std::size_t hash = init;
    for (std::size_t i = 0; i < key_size; ++i) {
        hash ^= key[i];
        hash *= multiplier;
    }

    return hash % BUCKET_INIT_COUNT;
}

template <class CharT>
bool keyEqual(const CharT* key_lhs, std::size_t key_size_lhs,
              const CharT* key_rhs, std::size_t key_size_rhs) {
    if (key_size_lhs != key_size_rhs) {
        return false;
    } else {
        return std::memcmp(key_lhs, key_rhs, key_size_lhs * sizeof(CharT)) == 0;
    }
}
}  // namespace hashRelative
}  // namespace myTrie

template <class CharT = char, class T = int, class KeySizeT = std::uint16_t>
class array_bucket {
   public:
    CharT* arr_buffer;
    size_t entry_count;
    size_t buffer_size;
    static const KeySizeT END_OF_BUCKET = std::numeric_limits<KeySizeT>::max();

    array_bucket() {
        entry_count = 0;
        // inited array_bucket: |END_OF_BUCKET|
        arr_buffer = (CharT*)std::malloc(sizeof(END_OF_BUCKET));
        std::memcpy(arr_buffer, &END_OF_BUCKET, sizeof(END_OF_BUCKET));
        buffer_size = sizeof(END_OF_BUCKET);
        std::cout << "init oldsize: " << buffer_size << std::endl;
    }

    bool is_end_of_bucket(const CharT* buffer) {
        return read_key_size(buffer) == END_OF_BUCKET;
    }

    size_t read_key_size(const CharT* buffer) {
        KeySizeT key_size;
        std::memcpy(&key_size, buffer, sizeof(KeySizeT));

        return (size_t)key_size;
    }

    // if found, return (entry_pos, true)
    // if not found, return (inserting_pos, false)
    std::pair<size_t, bool> find_in_bucket(const CharT* key, size_t keysize) {
        CharT* buffer_ptr = arr_buffer;
        size_t pos = 0;
        while (!is_end_of_bucket(buffer_ptr)) {
            size_t length = read_key_size(buffer_ptr);
            CharT* cmp_buffer_ptr = buffer_ptr + sizeof(KeySizeT);
            if (myTrie::hashRelative::keyEqual(key, keysize, cmp_buffer_ptr,
                                               length)) {
                return std::pair<size_t, bool>(pos, true);
            }
            // move ptr to next header, skip keysize, string, value
            buffer_ptr = buffer_ptr + sizeof(KeySizeT) + length + sizeof(T);
            pos += sizeof(KeySizeT) + length + sizeof(T);
        }
        return std::pair<size_t, bool>(pos, false);
    }

    size_t cal_arrbuffer_newsize(size_t keysize) {
        // when add a new key we need extra |size|string|value|
        size_t new_buffer_size = buffer_size;
        new_buffer_size +=
            sizeof(KeySizeT) + keysize * sizeof(CharT) + sizeof(T);
        return new_buffer_size;
    }

    T& append(const CharT* target, size_t keysize, T& value) {
        std::pair<size_t, bool> found = find_in_bucket(target, keysize);

        if (found.second) {
            // found the target, return the value reference
            return get_val_ref(found.first);

        } else {
            std::cout << "not found, start realloc" << std::endl;
            std::cout << (void*)found.first << " & " << (void*)arr_buffer
                      << std::endl;
            // not found, the buffer is full, need realloc
            size_t new_size = cal_arrbuffer_newsize(keysize);
            std::cout << "oldsize: " << buffer_size << std::endl;
            std::cout << "newsize: " << new_size << std::endl;

            CharT* new_buffer = (CharT*)(std::realloc(arr_buffer, new_size));
            if (new_buffer == nullptr) {
                std::cout << "realloc failed!!!!!!!!!1" << std::endl;
                exit(-1);
            }
            arr_buffer = new_buffer;
            buffer_size = new_size;

            append_impl(target, keysize, arr_buffer + found.first, value);

            return get_val_ref(found.first);
        }
    }

    void append_impl(const CharT* target, size_t keysize,
                     CharT* buffer_append_pos, T& value) {
        // append key_size
        std::memcpy(buffer_append_pos, &keysize, sizeof(KeySizeT));
        buffer_append_pos += sizeof(KeySizeT);

        // append the string
        std::memcpy(buffer_append_pos, target, keysize * sizeof(CharT));
        buffer_append_pos += keysize;

        // append the value
        std::memcpy(buffer_append_pos, &value, sizeof(T));
        buffer_append_pos += sizeof(T);

        // append the whole buffer end
        const auto end_of_bucket = END_OF_BUCKET;
        std::memcpy(buffer_append_pos, &end_of_bucket, sizeof(end_of_bucket));
    }

    T& get_val_ref(size_t val_pos) {
        CharT* entry_start_pos = arr_buffer + val_pos;
        size_t length = read_key_size(entry_start_pos);
        CharT* valAddress = entry_start_pos + sizeof(KeySizeT) + length;
        return *((T*)valAddress);
    }

    // debug
    void read_array_bucket() {
        std::cout
            << "\n\n-----------------------------\ncall read_array_bucket\n";
        CharT* buf = arr_buffer;
        while (!is_end_of_bucket(buf)) {
            size_t length = read_key_size(buf);
            std::cout << *((KeySizeT*)buf) << "|";

            for (size_t i = 0; i != length; i++) {
                std::cout << *(buf + sizeof(KeySizeT) + i) << "|";
            }
            // move ptr to next header, skip keysize, string, value
            buf = buf + sizeof(KeySizeT) + length;
            cout << (T)*buf << "|" << endl;
            buf = buf + +sizeof(T);
        }
        std::cout << "END read: |" << (int)*(buf) << "| finish " << std::endl;
        std::cout << "-----------------------------\n\n";
    }

    std::vector<std::pair<std::string, T>> get_item_in_array_bucket() {
        std::cout << "\n\n-----------------------------\ncall "
                     "get_item_in_array_bucket\n";
        std::vector<std::pair<std::string, T>> res;
        CharT* buf = arr_buffer;
        while (!is_end_of_bucket(buf)) {
            size_t length = read_key_size(buf);

            char* temp = (char*)malloc(length);
            std::memcpy(temp, buf + sizeof(KeySizeT), length);
            std::string item(temp);
            free(temp);

            // move ptr to next header, skip keysize, string, value
            buf = buf + sizeof(KeySizeT) + length;
            T v = *((T*)buf);

            res.push_back(std::pair<std::string, T>(item, v));
            buf = buf + sizeof(T);
        }
        return res;
    }
};

int main() {
    array_bucket<char, int, std::uint16_t> ab;
    ab.read_array_bucket();
    // std::string test1 = "haha";
    // int val1 = 1;
    // std::pair<char*, bool> res1 =
    //     ab.find_in_bucket(test1.c_str(), test1.size());
    // if (res1.second == false) {
    //     ab.append(test1.c_str(), test1.size(), res1.first, val1);
    // } else {
    //     std::cout << "find" << std::endl;
    // }
    std::string test1 = "hah1";

    int val1 = 1;
    cout << "testing: data:" << test1 << " size:" << test1.size() << endl;
    ab.append(test1.data(), test1.size(), val1) = 6;

    ab.read_array_bucket();

    std::string test2 = "hah2";
    cout << "testing: data:" << test2 << " size:" << test2.size() << endl;
    ab.append(test2.data(), test2.size(), val1) = 8;
    ab.read_array_bucket();

    std::string test3 = "hah3";
    cout << "testing: data:" << test3 << " size:" << test3.size() << endl;
    ab.append(test3.data(), test3.size(), val1) = 10;
    ab.read_array_bucket();

    char test[] = {'a', 'b', 'c'};
    int val2 = 5;
    cout << "testing: data:" << test << " size:" << sizeof(test) << endl;

    ab.append(test, sizeof(test), val2);

    // char test[] = {'a', 'b', 'c'};
    // std::pair<char*, bool> res = ab.find_in_bucket(test, sizeof(test));
    // // find
    // if (res1.second == false) {
    //     ab.append(test1.c_str(), test1.size(), res1.first, val1);
    // } else {
    //     std::cout << "find" << std::endl;
    // }

    ab.read_array_bucket();

    std::vector<std::pair<std::string, int>> res =
        ab.get_item_in_array_bucket();

    for (auto it = res.begin(); it != res.end(); it++) {
        cout << it->first << "-" << it->second << endl;
    }
}