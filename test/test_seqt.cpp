
#include "seqt.hpp"

#include <iostream> 

int main(int, char**) {
    using std::cin;

    mem m;

    const char * in = "ablkjasdrlkewjoasdigjierieoir2834rudkjfalksdhgj8q2934hrauisdjhfkasdjhf98234roiadsjsfa iu932r asdf ";

    for(int i = 0; i < strlen(in); i++) {
        m.read(in[i]);
    }

    while(cin) {
        auto c = cin.get();
        m.read(c);
    }
}