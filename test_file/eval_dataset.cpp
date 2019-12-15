#include <fstream>
#include <string>
#include <iostream>

using namespace std;

void checkDataset(string testing_dataset) {
    std::fstream f(testing_dataset);
    string url = "";
    uint32_t v = 0;;

    uint64_t charnum = 0;
    uint64_t count = 0;
    while (f >> url >> v) {
        charnum += url.size();
        count++;
    }
    cout << testing_dataset << ": " << charnum / count << endl;
}

int main() {
    checkDataset("dataset/id_yago/str_normal");
    checkDataset("dataset/id_lubm_640/str_normal");
}