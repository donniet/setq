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

    void for_each_thread(std::function<void(thread &)>);
    void reinforce(seqt);
    void disregard(seqt);

    seqt ordered(symbol);
    seqt unordered(symbol);

    seqt nil() const;

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
    std::vector<thread> threads_;
    std::deque<ident> recycled_;
    ident next_id_;

    std::map<symbol, ident> atom_index_;
    // an index which maps a seq A, to a set containing
    // all the seqs which have a single symbol appended to A
    std::map<tuple<ident, symbol>, ident> ordered_next_index_; 

    seqt create_atom(symbol c);
    seqt create_ordered(seqt prev, symbol c);
    seqt create_ordered(symbol c);
    seqt create_unordered(seqt prev, symbol c);
    seqt create_unordered(symbol c);

    bool find_in_unordered(seqt, symbol);

    mem::seqt splay(seqt, symbol);

    int compare(seqt, symbol);

    seqt_data & data(ident);
};

struct mem::seqt {
    friend mem::seqt_less;
    friend mem;

    struct iterator;
    struct pointer;

    // iterator begin();
    // iterator end();

    // iterator set_find(seqt const & s);

    seqt append(symbol, bool ordered_if_atom = true) const;
    seqt remove(symbol) const;
    symbol repr() const;
    bool expects(symbol);
    bool is_ordered() const;
    bool is_unordered() const;
    bool is_atom() const;
    bool is_nil() const;

    bool operator==(seqt const &) const;
private:
    ident id_;
    mem * owner_;

    seqt() : id_(0), owner_(nullptr) { }
    seqt(mem * owner) : id_(0), owner_(owner) { }
    seqt(ident id) : id_(id), owner_(nullptr) { }
    seqt(ident id, mem * owner) : id_(0), owner_(owner) { }
};

struct mem::seqt_data {
    ordinality ordinality_;
    ident left_;
    ident right_;
    symbol repr_;
};

struct mem::thread {
    bool advance(symbol);

    seqt & refers_to();
    seqt & visited(); // visited elements of the refers_to seqt
    seqt & unvisited(); // unvisited elements of the refers_to seqt
private:
    seqt visited_;
    seqt unvisited_;
};


mem::seqt mem::seqt::append(symbol c, bool ordered_if_atom) const {
    if(is_nil()) return owner_->create_atom(c);
    if(is_atom()) {
        if(ordered_if_atom) return owner_->create_ordered(*this, c);
        return owner_->create_unordered(*this, c);
    }
    if(is_ordered()) return owner_->create_ordered(*this, c);
    /* if(is_unordered()) */
    return owner_->create_unordered(*this, c);
}

bool mem::seqt::is_nil() const { return id_ == 0; }


bool mem::seqt::operator==(seqt const & a) const {
    return id_ == a.id_ && owner_ == a.owner_;
}

bool mem::seqt_less::operator()(mem::seqt const & a, mem::seqt const & b) const {
    if(a.owner_ < b.owner_) return true; // atoms are less than non-atoms
    if(b.owner_ < a.owner_) return false; 
    if(a.id_ < b.id_) return true;
    return false;
}

bool mem::thread::advance(symbol c) {
    if(!unvisited().expects(c))
        return false;

    // move this thread forward
    visited() = visited().append(c);
    unvisited() = unvisited().remove(c);
    return true;
}

mem::seqt mem::create_atom(symbol c) {
    return seqt(c);
}

mem::seqt mem::create_ordered(mem::seqt prev, symbol c) {
    auto a = create_atom(c);

    auto i = ordered_next_index_.find({prev.id_, c});
    if(i != ordered_next_index_.end()) 
        return seqt(i->second, this); // found it

    auto id = next_id_++;
    ordered_next_index_[{prev.id_, c}] = id;
    return seqt(id, this);
}

// returns an unordered seqt which is 
mem::seqt mem::splay(mem::seqt root, symbol c) {
    if(!root.is_unordered()) 
        return nil();

    mem::seqt ret = create_unordered(c);

    mem::seqt_data & dat = data(ret.id_);
    mem::seqt_data & rd = data(root.id_);
    ident * left = &dat.left_;
    ident * right = &dat.right_;



    while(root.is_unordered()) {
        if(c < rd.repr_) {
            *right = create_unordered(root.repr()).id_;
            mem::seqt_data & temp = data(*right);
            temp.right_ = rd.right_;
            rd = temp;
        } else if(c > root.repr()) {
            cr.left() = create_unordered(root.repr());
            cr = cr.left();
            cl->right() = root.right();
            root = *root.left();
        } else {
            break;
        }
    }

    return ret;
}

bool mem::find_in_unordered(mem::seqt in, symbol c) {
    auto a = create_atom(c);

    while(!in.is_nil() && in != a) {
        
    }
}

mem::seqt mem::create_unordered_without(mem::seqt s, symbol c) {
    if(!find_in_unordered(s, c)) // O(log N)
        return s;

    auto lt = create_unordered_lt(s, c); 
    auto gt = create_unordered_gt(s, c);

    return unordered_join(lt, gt);
}

/*

unordered seqt:
* find in O(log N)
* compare to other seqts in O(1) - O(log N + log M)


*/

mem::seqt mem::create_unordered(mem::seqt prev, symbol c) {
    auto a = create_atom(c);

    // is this atom already in prev?
    if(find_in_unordered(prev, a))
        return prev;

    // unorderd sets are binary trees, which are really ordered pairs
    // of binary trees. If we have an unordered set equal to {prev U {c}}
    // and prev already exists, then L = {x | x \in prev x < c} should exist
    auto lt = create_unordered_lt(prev, c);
    auto gt = create_unordered_gt(prev, c);

    return unordered_join(unordered_join(lt, c), gt);
}

void mem::for_each_thread(std::function<void(thread&)> func) {
    for(auto i = threads_.begin(); i != threads_.end(); i++)
        func(*i);
}

void mem::read(uint32_t in) {
    auto c = symbol(in);
    
    for_each_thread([c, this](thread & t) {
        if(t.advance(c)) {
            reinforce(t.visited());
        } else {
            disregard(t.visited());
            reinforce(t.visited().append(c));
        }
    });
    reinforce(ordered(c));
    reinforce(unordered(c));
}

int main(int, char**) {
    using std::cin;

    mem m;

    while(cin) {
        auto c = cin.get();
        m.read(c);
    }
}