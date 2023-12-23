#ifndef __SEQT_HPP__
#define __SEQT_HPP__

#include <numbers>
#include <functional>
#include <iostream>
#include <queue>
#include <map>
#include <tuple>

using std::tuple;
using std::tie;

#define MAX_UNICODE 0x10FFFF

struct symbol {
    uint32_t character_code_;

    operator uint32_t() const { return character_code_; }
    bool operator==(symbol const & r) const { return character_code_ == r.character_code_; }
    bool operator==(symbol && r) const { return character_code_ == r.character_code_; }
    bool operator!=(symbol const & r) const { return !(*this == r); }
    bool operator!=(symbol && r) const { return !(*this == r); }
    bool operator<(symbol const & r) const { return character_code_ < r.character_code_; }
    bool operator<(symbol && r) const { return character_code_ < r.character_code_; }

    bool is_unicode() const { return character_code_ <= MAX_UNICODE; }

    // special symbols
    static symbol nil() { return symbol{0x00FFFFFF}; }
    static symbol ordered() { return symbol{0x01FFFFFF}; }
    static symbol unordered() { return symbol{0x2FFFFFF}; }

    symbol(uint32_t c) : character_code_(c) { }
};

class mem {
    struct seqt;
    struct thread;

    typedef uint64_t ident;

    friend struct seqt;

    
    void process_threads(std::function<void(thread &)>);
    void reinforce(thread);
    void disregard(thread);

    seqt nil();

public:
    void read(uint32_t c);

private:
    typedef enum {
        ordered = 0,
        unordered = 1
    } ordinality;

    struct seqt_data;

    struct seqt_less {
        bool operator()(seqt const & a, seqt const & b) const;
    };

    // uses the raw operations
    seqt create_ordered(seqt prev, symbol c);
    seqt create_ordered(symbol c);
    seqt create_unordered(seqt prev, symbol c);
    seqt create_unordered(symbol c);

    bool find_in_unordered(seqt, symbol);

    seqt maximum_unordered(seqt);
    ident merge_unordered_ids(ident, ident);
    seqt splay_insert(seqt, symbol);
    seqt splay_remove(seqt, symbol);

    // raw operations
    seqt_data * get_data(ident);
    tuple<ident,seqt_data*> new_seqt(symbol, ordinality);

    // data members

    typedef std::vector<thread> thread_container;
    thread_container threads_;

    std::deque<ident> recycled_;
    std::vector<seqt_data> data_;
    ident next_id_;

    // indexes
    std::map<symbol, ident> unordered_atom_index_;
    std::map<symbol, ident> ordered_atom_index_;
    // std::map<ident, symbol> symbol_index_;

    std::map<tuple<ident,symbol>,ident> ordered_next_index_;

public:
    mem();
};

struct mem::seqt {
    friend mem::seqt_less;
    friend mem;

    struct iterator;
    struct pointer;

    // iterator begin();
    // iterator end();

    // iterator set_find(seqt const & s);

    seqt append(symbol, bool default_ordered = true) const;
    seqt remove(symbol) const;
    symbol repr() const;
    seqt expects(symbol);
    bool is_ordered() const;
    bool is_unordered() const;
    bool is_nil() const;

    operator bool() const { return !is_nil(); }
    bool operator==(seqt const &) const;
private:
    ident id_;
    mem * owner_;

    seqt() : id_(0), owner_(nullptr) { }
    seqt(mem * owner) : id_(0), owner_(owner) { }
    seqt(ident id) : id_(id), owner_(nullptr) { }
    seqt(ident id, mem * owner) : id_(id), owner_(owner) { }
};

struct mem::seqt_data {
    ordinality ordinality_;
    symbol repr_;
    ident left_;
    ident right_;
};

struct mem::thread {
    bool advance(symbol);

    seqt refers_to();
    seqt visited(); // visited elements of the refers_to seqt
    seqt unvisited(); // unvisited elements of the refers_to seqt
private:
    seqt visited_;
    seqt unvisited_;
public:
    thread(mem::seqt un) : visited_(0, un.owner_), unvisited_(un)  { }
    thread(mem::seqt vis, mem::seqt un) : visited_(vis), unvisited_(un) { } 
};

#endif // __SEQT_HPP__