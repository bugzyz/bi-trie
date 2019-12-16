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

int main() {
    std::cout << "starto!\n";

    string output_file = "slice_data/result";
    std::fstream result_output_f(output_file, std::ios::out | std::ios::app);

#ifdef TEST_YAGO
    string testing_dataset = "dataset/id_yago/cut_str_normal";
#else
    string testing_dataset = "dataset/id_lubm_640/str_normal";
#endif

    cout << "testing file: " << testing_dataset << endl;

    //------------testing--------------
    cout << "\n========unordered map result==========\n";

    std::string url;
    uint32_t v;
    uint64_t staUm;
    uint64_t endUm;
    uint64_t staTm;
    uint64_t endTm;

    std::fstream f(testing_dataset);
    boost::unordered_map<string, uint32_t> m1;
    boost::unordered_map<uint32_t, string> m2;

    staUm = get_time();

    double virt = 0.0;
    double res = 0.0;
    
    uint32_t count = 0;
    // uint64_t orignal_mem_start = get_start_cur_memory_clear();
    typedef struct mallinfo malloc_info;
    malloc_info mem_start, mem_end;

    mem_start = mallinfo();
    while (f >> url >> v) {
        m1[url] = v;
        m2[v] = url;
        count++;
    }
    mem_end = mallinfo();
    // uint64_t orignal_mem_end = get_end_cur_memory_clear();
    // cout << "original mem: " << orignal_mem_end - orignal_mem_start << " mb" << endl;
    cout << "original mem: " << mem_end - mem_start << " mb" << endl;
    
    endUm = get_time();

    f.close();

    cout << "unordered_map use time: usec: " << endUm - staUm << std::endl;

    vector<size_t> associativitys;
    vector<size_t> buckets;

    // analyse by elem_per_bucket
    // associativity should be 4, 8 for alignment
    associativitys.push_back(4);
    // associativitys.push_back(8);

    // analyse by bucket_num
    // buckets.push_back(463);
    // buckets.push_back(271);
    // buckets.push_back(101);
    // buckets.push_back(71);
    buckets.push_back(59);
    // buckets.push_back(31);
    // buckets.push_back(11);
    // buckets.push_back(5);

    // pair<Bucket, associativity>
    vector<pair<size_t, size_t>> configs;
    for (auto ass : associativitys) {
        for (auto buc : buckets) {
            configs.push_back(std::pair<size_t, size_t>(buc, ass));
        }
    }

    for (auto it = configs.begin(); it != configs.end(); it++) {
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

        size_t buc = it->first;
        size_t ass = it->second;
        cout << "testing: (" << buc << ", " << ass << ")" << endl;

        bi_trie<char, uint32_t> *hm_ptr = new bi_trie<char, uint32_t>();
        bi_trie<char, uint32_t> &hm = *hm_ptr;
        std::fstream f1(testing_dataset);

        staTm = get_time();
        int ind = 0;

        mem_start = mallinfo();
        while (f1 >> url >> v) {
            ind++;

            hm.insert_kv(url, v);
            if ((ind % 100000) == 0) cout << ind << endl;
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

        /*-----------------access performance calculating-----------------*/
        int64_t hm_k_total_time = 0;
        int64_t um_k_total_time = 0;

        // load key and value into independent vectors
        vector<string> string_vec;
        vector<uint32_t> id_vec;

        for(auto it = m1.begin();it!=m1.end();it++){
            string_vec.push_back(it->first);
            id_vec.push_back(it->second);
        }

        std::string url0;
        std::string url1;
        std::string url2;
        std::string url3;
        std::string url4;
        std::string url5;
        std::string url6;
        std::string url7;
        std::string url8;
        std::string url9;

        // checking:
        for (int i = 0; i < string_vec.size(); i = i + 10) {
            url0 = string_vec[i++];
            url1 = string_vec[i++];
            url2 = string_vec[i++];
            url3 = string_vec[i++];
            url4 = string_vec[i++];
            url5 = string_vec[i++];
            url6 = string_vec[i++];
            url7 = string_vec[i++];
            url8 = string_vec[i++];
            url9 = string_vec[i++];

            if(!(i < string_vec.size()))
                break;

            uint32_t gotfromhm;
            int64_t hm_get_start = get_time();
            gotfromhm = hm[url0];
            gotfromhm = hm[url1];
            gotfromhm = hm[url2];
            gotfromhm = hm[url3];
            gotfromhm = hm[url4];
            gotfromhm = hm[url5];
            gotfromhm = hm[url6];
            gotfromhm = hm[url7];
            gotfromhm = hm[url8];
            gotfromhm = hm[url9];
            int64_t hm_get_end = get_time();

            uint32_t gotfromum;
            int64_t um_get_start = get_time();
            gotfromum = m1[url0];
            gotfromum = m1[url1];
            gotfromum = m1[url2];
            gotfromum = m1[url3];
            gotfromum = m1[url4];
            gotfromum = m1[url5];
            gotfromum = m1[url6];
            gotfromum = m1[url7];
            gotfromum = m1[url8];
            gotfromum = m1[url9];
            int64_t um_get_end = get_time();

            int64_t um_used_time = um_get_end - um_get_start;
            int64_t hm_used_time = hm_get_end - hm_get_start;

            hm_k_total_time += hm_used_time;
            um_k_total_time += um_used_time;

            double cur_percent_k =
                (double)hm_used_time / (double)um_used_time * (double)100;

            if (max_percent_k < cur_percent_k) {
                max_percent_k = cur_percent_k;
            }
            if (min_percent_k > cur_percent_k) {
                min_percent_k = cur_percent_k;
            }
        }
        search_key_avg_time = (double)hm_k_total_time /  (double)count;
        cout << "search_key_avg_time: "
             <<  search_key_avg_time << endl;

        search_key_compared_with_baseline =
            (double)hm_k_total_time / (double)um_k_total_time * 100.0;
        cout << "search_key_compared_with_baseline: "
             << search_key_compared_with_baseline << "%" << endl;

        // checking:
        int64_t hm_v_total_time = 0;
        int64_t um_v_total_time = 0;

        uint32_t v0;
        uint32_t v1;
        uint32_t v2;
        uint32_t v3;
        uint32_t v4;
        uint32_t v5;
        uint32_t v6;
        uint32_t v7;
        uint32_t v8;
        uint32_t v9;

        for (int i = 0; i < id_vec.size(); i = i + 10) {
            v0 = id_vec[i++];
            v1 = id_vec[i++];
            v2 = id_vec[i++];
            v3 = id_vec[i++];
            v4 = id_vec[i++];
            v5 = id_vec[i++];
            v6 = id_vec[i++];
            v7 = id_vec[i++];
            v8 = id_vec[i++];
            v9 = id_vec[i++];

            if(!(i < id_vec.size()))
                break;

            std::string gotfromhm;
            int64_t hm_get_start = get_time();
            gotfromhm = hm[v0];
            gotfromhm = hm[v1];
            gotfromhm = hm[v2];
            gotfromhm = hm[v3];
            gotfromhm = hm[v4];
            gotfromhm = hm[v5];
            gotfromhm = hm[v6];
            gotfromhm = hm[v7];
            gotfromhm = hm[v8];
            gotfromhm = hm[v9];
            int64_t hm_get_end = get_time();

            std::string gotfromum;
            int64_t um_get_start = get_time();
            gotfromhm = m2[v0];
            gotfromhm = m2[v1];
            gotfromhm = m2[v2];
            gotfromhm = m2[v3];
            gotfromhm = m2[v4];
            gotfromhm = m2[v5];
            gotfromhm = m2[v6];
            gotfromhm = m2[v7];
            gotfromhm = m2[v8];
            gotfromhm = m2[v9];
            int64_t um_get_end = get_time();


            int64_t um_used_time = um_get_end - um_get_start;
            int64_t hm_used_time = hm_get_end - hm_get_start;

            hm_v_total_time += hm_used_time;
            um_v_total_time += um_used_time;

            double cur_percent_v =
                (double)hm_used_time / (double)um_used_time * (double)100;

            if (max_percent_v < cur_percent_v) {
                max_percent_v = cur_percent_v;
            }
            if (min_percent_v > cur_percent_v) {
                min_percent_v = cur_percent_v;
            }
        }

        search_value_avg_time = (double)hm_v_total_time / (double)count;
        cout << "search_value_avg_time: " << search_value_avg_time << endl;

        search_value_compared_with_baseline =
            (double)hm_v_total_time / (double)um_v_total_time * 100.0;
        cout << "search_value_compared_with_baseline: "
             << search_value_compared_with_baseline << "%" << endl;

        // correctness check
        if (test_and_print_wrong_test) {
            vector<string> wrong_search_key;
            vector<uint32_t> wrong_search_value;

            for (auto it = m1.begin(); it != m1.end(); it++) {
                if (it->second != hm[it->first]) {

                    string se = it->first;
                    uint32_t sey = it->second;
                    
                    uint32_t m1_res = m1[se];
                    uint32_t hm_res = hm[se];

                    cout << "wrong key search: " << it->first << endl;
                    cout << "get: " << hm[it->first] << endl;
                    wrong_search_key.push_back(it->first);
                }
            }

            // cout << "test key finish!\n";
            // cout << "wrong key_searching num: " << wrong_search_key.size()
            //      << endl;

            for (auto it = m2.begin(); it != m2.end(); it++) {
                if (it->second != hm[it->first]) {
                    cout << "wrong value search: " << it->first << endl;
                    cout << "get: " << hm[it->first] << endl;
                    wrong_search_value.push_back(it->first);
                }
            }

            // cout << "test value finish!\n";
            // cout << "wrong value_searching num: " << wrong_search_value.size()
            //      << endl;

            if (wrong_search_value.size() != 0 ||
                wrong_search_key.size() != 0) {
                cout << "testing failed\n";
                exit(0);
            }
        }

        if (manually_test) {
            while (cin >> url) {
                cout << "get value: " << hm[url] << endl;
            }
        }

        // Evaluation print out
        result_output_f << buc;
        result_output_f << ",";
        result_output_f << ass;
        result_output_f << ",";

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
    }
    result_output_f.close();
}