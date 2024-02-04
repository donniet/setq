#ifndef __UTF8_READER_HPP__
#define __UTF8_READER_HPP__

#include <iostream>

using std::istream;

struct utf8_reader {
    istream & is_;

    utf8_reader(istream & is) : is_(is) { }

    char32_t get() {
        static const unsigned char first_bits[] = {0xff, 0x7f, 0x3f, 0x1f, 0x0f};

        char c = is_.get();
        
        int bytes = 0;
        for(char d = c; d & 0x80; bytes++, d <<= 1) 
            {}

        char32_t ret = c & first_bits[bytes];
        for(; bytes > 0 && !is_.eof(); bytes--) {
            c = is_.get();
            ret <<= 7;
            ret |= c & 0x7f;
        }

        return ret;
    }
    bool eof() {
        return is_.eof();
    }

};


#endif