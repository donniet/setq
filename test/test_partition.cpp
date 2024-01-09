#include "partition.hpp"

#include <iostream>
#include <iomanip>

using namespace std;

int main(int ac, char ** av) {
    cout << setw(5) << " ";
    for(unsigned long n = 0; n < 10; n++) {
        cout << setw(5) << n;
    }
    cout << "\n";
    for(unsigned long k = 0; k < 10; k++) {
        cout << setw(5) << k;
        for(unsigned long n = 0; n < 10; n++) {
            cout << setw(5) << partition_identical_n_distinct_k(n, k);
        }
        cout << "\n";
    }
    cout << endl;

    // cout << "P(4,4): " << partition_identical_n_distinct_k(4, 5) << endl;

    return 0;
}



