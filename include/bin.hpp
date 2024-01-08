

#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <stdint.h>

using namespace std;

/*
    binary version of sequence-set
 */

struct transition {
    size_t from;
    size_t to;
};

struct binseqt {
    vector<size_t> weights;
    size_t count;

    void read(bool next_bit);

    template<typename Word>
    void read_word(Word w) {
        // big-endian
        Word mask = 0x1 << (sizeof(Word) * 8 - 1);
        for(; mask > 0; mask >>= 1) {
            read(mask & w);
        }
    }

    binseqt() : count(0), weights(2) {
        weights[0] = 0;
        weights[1] = 0;
    }
};

void binseqt::read(bool bit) {
    // increment the bit weight
    
}