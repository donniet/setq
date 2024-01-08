#include "test.hpp"
#include "manual.hpp"

#include <iostream>
#include <string>

using std::cin;
using std::cout;
using std::string;


int main(int ac, char * av[]) {
    seqt s;

    if(ac > 1) { 
        while(cin) {
            uint32_t c = (uint32_t)cin.get();
            if(c < 0x10000)
                s.read(c);
        }
    } else {
        // test mode
        string test = "hello!  this is a test of a lot of text to see what happens and if we can detect any patterns in the output.";

        for(char c : test) {
            s.read(c);
        }
    }

    // s.cleanup(3);

    s.dump(cout);
    cout.flush();
    
    return 0;
}