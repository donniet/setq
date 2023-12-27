#include "test.hpp"
#include "manual.hpp"

#include <iostream>

using std::cin;
using std::cout;

int main(int ac, char * av[]) {
    seqt s;

    while(cin) {
        s.read(cin.get());
    }

    s.dump(cout);
    cout.flush();
    
    return 0;
}