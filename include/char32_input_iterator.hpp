#ifndef __CHAR32_INPUT_ITERATOR__
#define __CHAR32_INPUT_ITERATOR__


#include <fstream>
#include <tuple>

using namespace std;

struct char32_input_iterator {
    istream * is_;
    char * buf_;
    size_t b_;
    size_t e_;
    enum encoding_type {
        binary = 0, ascii, utf8 //, utf16, utf32
    } encoding_;
    size_t bytes_read_;
    char32_t current_utf8_char;
    size_t symbols_read_;

    static const size_t buffer_size_ = 2*BUFSIZ;
    static const char32_t maximum_utf8 = 0x80000000;

    size_t fill_buffer() {
        size_t maximum_bytes = get_remaining_buffer();
        bytes_read_ = 0;

        if(maximum_bytes == 0) 
            return 0;
        
        tuple<size_t,size_t> first_read =  {e_, e_ + maximum_bytes};
        tuple<size_t,size_t> second_read = { 0, 0 };
        if(get<1>(first_read) > buffer_size_) {
            get<1>(first_read) = buffer_size_;
            get<1>(second_read) = e_ + maximum_bytes - buffer_size_;
        }
      
        is_->read(&buf_[e_], get<1>(first_read) - get<0>(first_read));
        bytes_read_ += is_->gcount(); 
        e_ = (e_ + bytes_read_) % buffer_size_;

        if(is_->eof() || get<1>(second_read) == get<0>(second_read))  // || bytes_read_ < get<1>(first_read) - get<0>(first_read)
            return bytes_read_;

        is_->read(&buf_[0], get<1>(second_read) - get<0>(second_read));
        bytes_read_ += is_->gcount();
        e_ = is_->gcount() % buffer_size_;

        return bytes_read_;
    }

    size_t get_remaining_buffer() const {
        if (b_ == e_) return buffer_size_;
        if (b_ <  e_) return e_ - b_;
        else          return e_ + buffer_size_ - b_ - 1;
    }

    encoding_type process_start_bytes() {
        static const unsigned char utf8_bom[] = { 0xef, 0xbb, 0xbf };
        
        for(int i = 0; i < sizeof(utf8_bom)/sizeof(utf8_bom[0]); i++) {
            if(current_byte() != utf8_bom[i])
                return ascii;
            
            move_forward(); // can ignore the return here (i think)
        }
        return utf8;
    }

    char current_byte() const {
        return buf_[b_];
    }

    bool move_forward() {
        b_ = (b_ + 1) % buffer_size_;
        if(b_ == e_) {
            if(fill_buffer() == 0)
                return false;
        }

        return true;
    }

    // assumes encoding_ == utf8
    char32_t read_utf8() {
        static const unsigned char utf8_b1_masks[] = {0x7f, 0x1f, 0x0f, 0x09, 0x03, 0x01};
        static const unsigned char utf8_b1_bits[]  = {0x00, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc};
        static const int  shift_by[]      = {   0,    5,    4,    3,    2,    1};
        static const unsigned char trailing_bit_mask = 0x3f;
        static const unsigned char trailing_bit_head = 0x80;
        static int        trailing_bits = 6;
        static const int utf8_max_bytes = sizeof(utf8_b1_masks)/sizeof(utf8_b1_masks[0]);

        if(current_utf8_char < maximum_utf8)
            return current_utf8_char;

        char32_t ret = 0;
        
        size_t starting_b = b_;
        int bytes_read = 0;
        char c = current_byte();
        int toread = 0;

        for(; toread < utf8_max_bytes; toread++) {
            // do the start bits signify toread bytes remaining?
            if((c & (~utf8_b1_masks[toread])) == utf8_b1_bits[toread])
                break;
        }

        if(toread >= utf8_max_bytes)
            goto invalid_file_format;

        ret |= (utf8_b1_bits[toread+1] & c);

        for(; bytes_read < toread && move_forward(); bytes_read++) {
            c = current_byte();
            if((c & (~trailing_bit_mask)) != trailing_bit_head) 
                break;

            ret <<= trailing_bits; // 6 bits 
            ret |= (trailing_bit_mask & c);
        }

        if(bytes_read < toread)
            goto invalid_file_format;
    
        current_utf8_char = ret;
        symbols_read_++;
        return ret;

    invalid_file_format:
        current_utf8_char = maximum_utf8;
        b_ = starting_b; // reset our start cursor
        encoding_ = binary;
        return (char32_t)current_byte();
    }

    char32_t operator*() {
        switch(encoding_) {
        case utf8:
            return read_utf8();
        case binary:
        case ascii:
            return (char32_t)current_byte();
        }

        return 0;
    }

    char32_input_iterator & operator++() {
        switch(encoding_) {
        case utf8:
            if(current_utf8_char >= maximum_utf8) {
                // this can change the encoding_
                read_utf8();
            }
            // it could be that this file is not utf8 and the read_utf8 function errored
            // and changed our encoding.  If that's the case we are now encoding_ == binary
            if(encoding_ == utf8) {
                // mark the current character as read
                // so that the next read call will actually read
                current_utf8_char = maximum_utf8;
                break;
            }
            // otherwise fallthrough
        case binary:
        case ascii:
            move_forward();
            break;
        }
        return *this;
    }

    bool eof() const {
        if(is_ == nullptr) return true;
        if(!is_->eof()) return false;

        return buffer_size_ - get_remaining_buffer() == 0;
    }

    bool operator==(char32_input_iterator const & r) const {
        if(eof() && r.eof()) return true;
        if(eof() || r.eof()) return false;
        
        
        return is_ == r.is_ 
            &&    is_->tellg() + (long long)(buffer_size_ -   get_remaining_buffer()) 
            ==  r.is_->tellg() + (long long)(buffer_size_ - r.get_remaining_buffer());
    }

    inline bool operator!=(char32_input_iterator const & r) const {
        return !(*this == r);
    }

    char32_input_iterator(istream & is, encoding_type encoding = binary) 
        : is_(&is), b_(0), e_(0), 
          current_utf8_char(maximum_utf8), symbols_read_(0),
          encoding_(encoding)
    { 
        buf_ = new char[buffer_size_];
        fill_buffer();
        if(encoding == utf8) {
            process_start_bytes();
        }
        // TODO: should we override the encoding type or have an autodetect?
    }
    char32_input_iterator() 
        : is_(nullptr), b_(0), e_(0), current_utf8_char(maximum_utf8), symbols_read_(0), encoding_(binary), buf_(nullptr)
    { }

    char32_input_iterator(char32_input_iterator const &) = delete;
    char32_input_iterator(char32_input_iterator && r) 
        : is_(r.is_), b_(r.b_), e_(r.e_), 
          current_utf8_char(r.current_utf8_char), symbols_read_(r.symbols_read_),
          encoding_(r.encoding_), buf_(r.buf_)
    {
        r.b_ = 0;
        r.e_ = 0;
        r.current_utf8_char = 0;
        r.symbols_read_ = 0;
        r.encoding_ = binary;
        r.buf_ = nullptr;
    }
    ~char32_input_iterator() {
        delete [] buf_;
    }
};


#endif