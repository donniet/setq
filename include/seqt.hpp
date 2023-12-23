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

class seqt {
public:

    struct node;
    struct thread;

    typedef uint64_t ident;

    friend struct node;

    
    void process_threads(std::function<void(thread &)>);
    void reinforce(thread);
    void disregard(thread);

    node nil();

public:
    void read(uint32_t c);

public:
    typedef enum {
        ordered = 0,
        unordered = 1
    } ordinality;

    struct seqt_data;

    struct seqt_less {
        bool operator()(node const & a, node const & b) const;
    };

    // uses the raw operations
    node create_ordered(node prev, symbol c);
    node create_ordered(symbol c);
    node create_unordered(node prev, symbol c);
    node create_unordered(symbol c);

    bool find_in_unordered(node, symbol);

    node maximum_unordered(node);
    ident merge_unordered_ids(ident, ident);
    node splay_insert(node, symbol);
    node splay_remove(node, symbol);

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
    seqt();
};

struct seqt::node {
    friend seqt::seqt_less;
    friend seqt;

    struct iterator;
    struct pointer;

    // iterator begin();
    // iterator end();

    // iterator set_find(seqt const & s);

    node append(symbol, bool default_ordered = true) const;
    node remove(symbol) const;
    symbol repr() const;
    node expects(symbol);
    bool is_ordered() const;
    bool is_unordered() const;
    bool is_nil() const;

    operator bool() const { return !is_nil(); }
    bool operator==(node const &) const;
public:
    ident id_;
    seqt * owner_;

    node() : id_(0), owner_(nullptr) { }
    node(seqt * owner) : id_(0), owner_(owner) { }
    node(ident id) : id_(id), owner_(nullptr) { }
    node(ident id, seqt * owner) : id_(id), owner_(owner) { }
};

struct seqt::seqt_data {
    ordinality ordinality_;
    symbol repr_;
    ident left_;
    ident right_;
};

struct seqt::thread {
    bool advance(symbol);

    node refers_to();
    node visited(); // visited elements of the refers_to seqt
    node unvisited(); // unvisited elements of the refers_to seqt
public:
    node visited_;
    node unvisited_;
public:
    thread(seqt::node un) : visited_(0, un.owner_), unvisited_(un)  { }
    thread(seqt::node vis, seqt::node un) : visited_(vis), unvisited_(un) { } 
};

#endif // __SEQT_HPP__