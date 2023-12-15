#include <iostream>
#include <functional>
#include <vector>
#include <map>

typedef uint32_t symbol;
typedef uint64_t identifier;

struct mem;
struct thread;

struct seqt {
    friend struct mem;

    seqt();
    bool operator==(seqt) const;
    operator bool() const;

private:
    seqt(identifier id);

    identifier id_; 
};

// the null seqt
seqt::seqt() : id_(0) { }
seqt::seqt(identifier id) : id_(id) { } 
seqt::operator bool() const { return id_ == 0; }
bool seqt::operator==(seqt s) const { return id_ == s.id_; }


// represents a "position" in a seqt by tracking the subsequences before and after the current symbol
// ______*          head() sequence so far (* marks the predicted symbol, not in head) 
//       _______    tail()
// a b c d e f g    following_
//       ^
struct thread {
    friend class mem;

    thread(seqt following);

    seqt head() const;  // the visited elements of this seqt
    seqt tail() const;  // the unvisited elements of this seqt
    operator seqt() const { return following_; }
private:
    seqt tail_;
    seqt following_;
};

thread::thread(seqt following)
    : following_(following), tail_(following)
{ }

seqt thread::tail() const { return tail_; }



struct mem {
    void read(symbol);

    void reinforce(seqt);
    void disregard(seqt);

    seqt guess();

private:
    size_t weight(thread);
    void follow(seqt);
    bool is_following(seqt);

private:
    // storage and creation methods
    void for_each_thread(std::function<void(thread &)>);
    
    seqt atom(symbol c);
    seqt create_set(seqt a, seqt b);
    seqt create_seq(seqt a, seqt b);
    void increment_weight(seqt);
    bool advance_thread(thread& t, symbol c);
    identifier get_atom_id(symbol c);

    typedef enum {
        nil = 0, atom, sequence, set
    } seqt_type;

    struct seqt_data {
        seqt_type type;
        identifier id;
        //TODO: look into splay trees as a way to efficiently remove an item?
        identifier left;    // used as a binary tree if this is a set, 0 if atom
        identifier right;   // used as a linked list if this is a seq, 0 if atom
        identifier first;   // this is a reference to the seqt that is first (either top of the binary tree, or first in the sequence)
                            // this MUST NOT be nil if this is a sequence or set
        symbol c;           // set to 0 if this isn't an atom
        size_t weight;
    };

    struct tree_node {
        tree_node * left;
        tree_node * right;
        identifier data;
    };

    std::vector<seqt_data> data_;
    std::map<symbol, identifier> atom_index_;
    identifier next_id_;
public:
    mem();
};

mem::mem() 
    : next_id_(1), data_(1)
{
    data_[0] = seqt_data{nil,0,0,0,0,0,0};
    atom_index_[0] = 0;
}

bool mem::advance_thread(thread& t, symbol c) {
    if(c == 0) return false;

    identifier atom_id = get_atom_id(c);

    seqt_data * tail_data = &data_[t.tail_.id_];
    seqt_data * first_data;

    // check the symbol of tail_data
    switch(tail_data->type) {
    case nil:
        return false;
    case atom:
        if(tail_data->id == atom_id) {
            t.tail_ = seqt();
            return true;
        }
        return false;
    case sequence:
        if(tail_data->first == atom_id) {
            t.tail_ = seqt(tail_data->right);
            return true;
        }
        // some kind of recursion here I think...

        return false;
    case set:
        // we have to look in the binary tree to find atom_id
        // tricky part is any node could be sequence or a set
        // but I think that since we are looking for a symbol
        // we can just look for atoms in this tree maybe?
        
        // I don't like this.  Maybe use a splay tree to bubble the found item
        // to the top.  Then make a new seqt entry for the splay tree with
        // the root node removed

        

        tree_node *  new_tail = nullptr;
        tree_node ** new_tail_inserter = &new_tail;

        first_data = tail_data;
        while(first_data->first != 0) {
            if(first_data->first == atom_id) {
                // found it!
                // somehow "advance" this set
                return true;
            }
            *new_tail_inserter = new tree_node{nullptr, nullptr, first_data->first};

            if(first_data->first < atom_id) {
                new_tail_inserter = &(*new_tail_inserter)->left;
            } else {

            }
        }
    }

    return false;
}

void mem::read(symbol c) {
    // check which tracked sequences expected this symbol
    // collapse the stack
    for_each_thread([c, this](auto & s) {
        if(advance_thread(s, c)) {
            // c advances the seqt
            if(!s.tail()) {
                // this is empty
                reinforce(s);
            }
        } else {
            // c doesn't advance the seqt
            // should we investigate whether s and c are part of the same class?
            follow(create_set(s.tail(), atom(c)));
            // and investigate if part of the head of this seqt and c make a sequence?
            follow(create_seq(s.head(), atom(c)));
        }
    });
    // and follow this sequence
    follow(atom(c));
}

// make a guess as to what comes next
seqt mem::guess() {
    
}

void mem::reinforce(seqt s) 
{
    increment_weight(s);
}

void mem::follow(seqt s) 
{
    // are we already following this seqt?
    if(is_following(s))
        return;

    
}


int main(int, char**) {
    std::cout << "Hello World!\n";
    return 0;
}