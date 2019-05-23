#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace std;

void process_mem_usage(string part_name) {
    using namespace std;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
        tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
        stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
        starttime >> vsize >> rss;  // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                        1024;  // in case x86-64 is configured to use 2MB pages

    std::cout << part_name << " : vm_usage: " << vsize / 1024 / 1024 << " res: " << rss * page_size_kb
              << std::endl;
}


uint64_t get_cur_memory_non_clear() {
    using namespace std;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
        tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
        stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
        starttime >> vsize >> rss;  // don't care about the rest

    stat_stream.close();
    return vsize;
}

#include <malloc.h>
#include <stdlib.h>
#include <string.h>

static void display_mallinfo() {
    struct mallinfo mi;
    mi = mallinfo();

    printf("Total non-mmapped bytes (arena):       %d\n", mi.arena);
    printf("# of free chunks (ordblks):            %d\n", mi.ordblks);
    printf("# of free fastbin blocks (smblks):     %d\n", mi.smblks);
    printf("# of mapped regions (hblks):           %d\n", mi.hblks);
    printf("Bytes in mapped regions (hblkhd):      %d\n", mi.hblkhd);
    printf("Max. total allocated space (usmblks):  %d\n", mi.usmblks);
    printf("Free bytes held in fastbins (fsmblks): %d\n", mi.fsmblks);
    printf("Total allocated space (uordblks):      %d\n", mi.uordblks);
    printf("Total free space (fordblks):           %d\n", mi.fordblks);
    printf("Topmost releasable block (keepcost):   %d\n", mi.keepcost);
}

int operator-(struct mallinfo after, struct mallinfo before){
    return ((after.uordblks - before.uordblks) + (after.hblkhd - before.hblkhd))/1024/1024;
}

const uint64_t calloc_page_size = 1024;
uint64_t get_start_cur_memory_clear() {
    uint64_t start_mem = get_cur_memory_non_clear();
    uint64_t idle_mem = 0;
    cout << "start mem: " << start_mem / 1024/ 1024 <<endl;
    int counter = 0;
    while(true) {
        char* idle_page = (char*)calloc(calloc_page_size, sizeof(char));
        idle_mem += calloc_page_size;
        counter++;
        if((get_cur_memory_non_clear()/1024/1024) > (start_mem/1024/1024)){
            break;
        }
    }
    cout << "end mem: " << get_cur_memory_non_clear()  / 1024/ 1024 <<endl;
    cout << "malloc time: " << counter << " return mem: " << (get_cur_memory_non_clear()) / 1024 / 1024 << endl;
    return (get_cur_memory_non_clear()) / 1024 / 1024;
}

uint64_t get_end_cur_memory_clear() {
    uint64_t start_mem = get_cur_memory_non_clear();
    uint64_t idle_mem = 0;
    cout << "start mem: " << start_mem / 1024/ 1024 <<endl;
    int counter = 0;
    while(true) {
        char* idle_page = (char*)calloc(calloc_page_size, sizeof(char));
        idle_mem += calloc_page_size;
        counter++;
        if((get_cur_memory_non_clear()/1024/1024) > (start_mem/1024/1024)){
            break;
        }
    }
    cout << "end mem: " << get_cur_memory_non_clear()  / 1024/ 1024 <<endl;
    cout << "malloc time: " << counter << " return mem: " << (get_cur_memory_non_clear() - idle_mem) / 1024 / 1024 << endl;
    return (get_cur_memory_non_clear() - idle_mem) / 1024 / 1024;
}
