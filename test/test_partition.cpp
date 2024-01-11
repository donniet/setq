#include "partition.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#define WIDTH 7

double binomial(double n, double k) {
    return tgamma(n+1) / tgamma(k+1) / tgamma(n-k+1);
}

double lbinomial(double n, double k) {
    return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

int main(int ac, char ** av) {
    cout << setw(WIDTH) << " ";
    for(unsigned long n = 0; n < 13; n++) {
        cout << setw(WIDTH) << n;
    }
    cout << "\n";
    for(unsigned long k = 1; k < 13; k++) {
        cout << setw(WIDTH) << k;
        for(unsigned long n = 0; n < 13; n++) {
            cout << setw(WIDTH) << exp(lbinomial(n+k-1, k-1)); // binomial
        }
        cout << "\n";
    }
    cout << endl;

    // cout << "P(4,4): " << partition_identical_n_distinct_k(4, 5) << endl;

    return 0;
}



