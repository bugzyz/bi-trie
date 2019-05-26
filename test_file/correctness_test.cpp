bool test_and_print_wrong_test = true;
bool manually_test = false;

#include "../unified_impl/test_switcher.hpp"
#include "../util/my_memory.hpp"

#include <boost/unordered_map.hpp>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>

#include <stdint.h>
#include <sys/time.h>

using namespace std;

int main(int argc, char *argv[]) {
    std::cout << "starto!\n";

    if (argc!=2) {
        cout << "enter test file please\n";
        exit(0);
    }
    string testing_dataset = string(argv[1]);

    string output_file = "correctness_result";
    std::fstream result_output_f(output_file, std::ios::out | std::ios::app);


    cout << "testing file: " << testing_dataset << endl;

    // Evaluation record
    double constructed_time = 0;

    double search_key_avg_time = 0;
    double search_value_avg_time = 0;

    unsigned int search_key_compared_with_baseline = 0;
    unsigned int search_value_compared_with_baseline = 0;

    double max_percent_k = 0.0;
    double min_percent_k = 100.0;

    double max_percent_v = 0.0;
    double min_percent_v = 100.0;

    int mem_used = 0;

    zyz_trie::bi_trie<char, uint32_t> *hm_ptr = new zyz_trie::bi_trie<char, uint32_t>(59, 8);
    zyz_trie::bi_trie<char, uint32_t> &hm = *hm_ptr;
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
    uint64_t staTm;
    uint64_t endTm;
    staTm = get_time();

    mem_start = mallinfo();
    int ind = 0;
    while (f1 >> url >> v) {
        hm.insert_kv(url, v);
        if ((++ind % 100000) == 0) cout << ind << endl;
    }
    // cleaning
    hm.storage_resize();
    mem_end = mallinfo();
    mem_used = mem_end - mem_start;
    cout << "opt mem: " << mem_used << " mb" << endl;
    endTm = get_time();
    constructed_time = (double)(endTm - staTm) / (double)1000 / (double)1000;

    f1.close();

    /*-----------------insert performance calculating-----------------*/
    cout << "\n==================\n";
    cout << "constructing time: " << constructed_time << endl;

    std::ifstream f2(testing_dataset, std::ios::in);
    if(f2) {
        cout << "open f1 successfully\n";
    } else {
        cout << "open f1 failed\n";
        exit(0);
    }
    int ind2 = 0;
    while (f2 >> url >> v) {
        if(hm[url] != v) {
            cout << "get wrong result of key: " << url << " answer:" << v << " get: " << hm[url] << endl;
            exit(0);
        }
        if(hm[v] != url) {
            cout << "get wrong result of value: " << v << " answer:" << url << " get: " << hm[v] <<  endl;
            exit(0);
        }
        if ((++ind2 % 100000) == 0) cout << ind2 << endl;
    }
    f2.close();

    cout << "test pass!" << endl;

    // Evaluation print out
    result_output_f << testing_dataset << ",";
#ifdef CUCKOO_THEN_EXPAND
    result_output_f << "true,";
#else
    result_output_f << "false,";
#endif
#ifdef ID_STR_NON_OPT
    result_output_f << "true,";
#else
    result_output_f << "false,";
#endif
#ifdef STR_ID_NON_OPT
    result_output_f << "true,";
#else
    result_output_f << "false,";
#endif
#ifdef BURST_NON_OPT
    result_output_f << "true,";
#else
    result_output_f << "false,";
#endif

    char *buffer = (char *)malloc(1024);
    int buffer_len = 0;
    string outside_memory_mes = "";
    string inside_memory_mes = "";
    string inside_time_mes = "";
    string outside_time_mes = "";

    outside_memory_mes = string(buffer, sprintf(buffer, "%d,", mem_used));

    inside_memory_mes = hm.memory_usage_report();

    outside_time_mes = string(buffer, sprintf(buffer, "%.2lf,%.2lf,%.2lf,%u,%u,%.2lf,%.2lf,%.2lf,%.2lf,",
                    constructed_time, search_key_avg_time,
                    search_value_avg_time, search_key_compared_with_baseline,
                    search_value_compared_with_baseline, max_percent_k,
                    min_percent_k, max_percent_v, min_percent_v));

    inside_time_mes = hm.time_cost_report();

    result_output_f << outside_memory_mes << inside_memory_mes
                    << outside_time_mes << inside_time_mes << endl;
    result_output_f.flush();
    free(buffer);
    result_output_f.close();
}