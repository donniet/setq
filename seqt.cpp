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

    // uses the raw operations
    seqt create_ordered(seqt prev, symbol c);
    seqt create_ordered(symbol c);
    seqt create_unordered(seqt prev, symbol c);
    seqt create_unordered(symbol c);

    bool find_in_unordered(seqt, symbol);

    mem::seqt splay(seqt, symbol);

    // raw operations
    seqt_data * get_data(ident);
    tuple<ident,seqt_data*> new_seqt(symbol, ordinality);

    // data members

    std::vector<thread> threads_;
    std::deque<ident> recycled_;
    std::vector<seqt_data> data_;
    ident next_id_;

    std::map<symbol, ident> atom_index_;
    std::map<ident, symbol> symbol_index_;

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
    symbol repr_;
    ident left_;
    ident right_;
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

bool mem::seqt::expects(symbol c) {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return false;

    if(data->ordinality_ == unordered) {
        // if we are unordered, then we expect something that we haven't already seen
        return !owner_->find_in_unordered(*this, c);
    }

    // if we are unordered, then we should see if c 
    // we need to look through all the seqt's that have c as a repr_
    // and see if we are the prev of any of them.
    // I think we should create an unordered seqt for each ordered seqt
    // that contains ordered pairs of adjacent repr_ maybe?
}

mem::seqt mem::seqt::append(symbol c, bool default_ordered) const {
    if(is_nil()) {
        if(default_ordered) return owner_->create_ordered(c);
        return owner_->create_unordered(c);
    }
    if(is_ordered()) return owner_->create_ordered(*this, c);
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

mem::seqt mem::create_ordered(mem::seqt prev, symbol c) {
    ident id;
    mem::seqt_data * temp;

    mem::seqt_data * p = get_data(prev.id_);
    tie(id, temp) = new_seqt(c, ordered);
    if(p == nullptr) {
        return seqt(id, this);
    }
    temp->left_ = prev.id_;
    return seqt(id, this);
}

// returns an unordered seqt which is 
mem::seqt mem::splay(mem::seqt root, symbol c) {
    mem::seqt_data * cur = get_data(root.id_);

    ident ret_id;
    mem::seqt_data * ret;
    tie(ret_id, ret) = new_seqt(c, unordered);
    if(cur == nullptr)
        return seqt(ret_id, this);

    // let's just create new nodes for now.  
    // maybe we can filter them later.

    ident * left_insert = &ret->left_;
    ident * right_insert = &ret->right_;

    ident temp_id = root.id_;
    mem::seqt_data * temp;

    while(cur != nullptr) {
        if(cur->ordinality_ != unordered) {
            // insert ordered seqts on the right
            *right_insert = temp_id;
            break;
        }
        
        if(cur->repr_ < c) {
            // we should move cur to the left
            tie(*left_insert, temp) = new_seqt(cur->repr_, unordered);
            temp->left_ = cur->left_;
            left_insert = &temp->right_;
            temp_id = cur->right_;
            cur = get_data(cur->right_);
        } 
        else if(c < cur->repr_) {
            tie(*right_insert, temp) = new_seqt(cur->repr_, unordered);
            temp->right_ = cur->right_;
            right_insert = &temp->left_;
            temp_id = cur->left_;
            cur = get_data(cur->left_);
        }
        else {
            *left_insert = cur->left_;
            *right_insert = cur->right_;
            break;
        }
    }

    return seqt(ret_id, this);
}

bool mem::find_in_unordered(mem::seqt root, symbol c) {
    mem::seqt_data * cur = get_data(root.id_);

    while(cur != nullptr && cur->ordinality_ == unordered) {
        if(cur->repr_ < c)
            cur = get_data(cur->right_);
        else if(c < cur->repr_)
            cur = get_data(cur->left_);
        else
            return true;
    }

    return false;
}

/*

unordered seqt:
* find in O(log N)
* compare to other seqts in O(1) - O(log N + log M)


*/

mem::seqt mem::create_unordered(mem::seqt prev, symbol c) {
    seqt s = splay(prev, c);
    if(s.is_nil()) {
        int id;
        seqt_data * _;
        tie(id, _) = new_seqt(c, unordered);
        return seqt(id, this);
    }
    return s;
}

mem::seqt mem::create_ordered(symbol c) {
    ident id;
    seqt_data * _;
    tie(id, _) = new_seqt(c, ordered);
    return seqt(id, this);
}

mem::seqt mem::create_unordered(symbol c) {
    ident id;
    seqt_data * _;
    tie(id, _) = new_seqt(c, unordered);
    return seqt(id, this);
}

tuple<mem::ident, mem::seqt_data*> mem::new_seqt(symbol c, ordinality ord) {
    ident id;
    seqt_data * data;
    if(!recycled_.empty()) {
        id = recycled_.back();
        recycled_.pop_back();
        data = get_data(id);

        data->ordinality_ = ord;
        data->left_ = 0;
        data->right_ = 0;
        data->repr_ = c;
    } else {
        id = next_id_++;
        data = &*data_.insert(data_.begin() + id, seqt_data{
            ord, c, 0, 0
        }); 
    }
    return {id, data};
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
    reinforce(create_ordered(c));
    reinforce(create_unordered(c));
}

int main(int, char**) {
    using std::cin;

    mem m;

    while(cin) {
        auto c = cin.get();
        m.read(c);
    }
}