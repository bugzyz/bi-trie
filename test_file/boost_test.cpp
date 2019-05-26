#include "../util/my_timer.hpp"
#include "../util/my_memory.hpp"

#include <boost/unordered_map.hpp>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>

#include <stdint.h>
#include <sys/time.h>

int main(int argc, char *argv[]) {
    std::cout << "starto!\n";

    if (argc!=2) {
        cout << "enter test file please\n";
        exit(0);
    }
    string testing_dataset = string(argv[1]);

    string output_file = "correctness_result";

    cout << "testing file: " << testing_dataset << endl;

    // Evaluation record
    double constructed_time = 0;

    int mem_used = 0;

    boost::unordered_map<string, uint32_t> m1;
    boost::unordered_map<uint32_t, string> m2;
    std::ifstream f1(testing_dataset);

    if(f1) {
        cout << "open f1 successfully\n";
    } else {
        cout << "open f1 failed\n";
        exit(0);
    }
    uint32_t count = 0;
    // uint64_t orignal_mem_start = get_start_cur_memory_clear();
    typedef struct mallinfo malloc_info;
    malloc_info mem_start, mem_end;

    std::string url;
    uint32_t v;
    uint64_t staUm;
    uint64_t endUm;
    staUm = get_time();

    mem_start = mallinfo();
    int ind = 0;
    while (f1 >> url >> v) {
        m1[url] = v;
        m2[v] = url;
        if ((++ind % 100000) == 0) cout << ind << endl;
    }
    mem_end = mallinfo();

    cout << "original mem: " << mem_end - mem_start << " mb" << endl;
    
    endUm = get_time();

    f1.close();

    cout << "unordered_map use time: usec: " << endUm - staUm << std::endl;
}