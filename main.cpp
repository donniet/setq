#include <iostream>
#include <functional>
#include <vector>
#include <map>
#include <list>
#include <set>
#include <exception>

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
        size_t weight;
        
        // atom parameters
        symbol c;

        // non-atom parameters
        identifier repr;

        // list parameters
        identifier seq_next;
        
        // tree parameters
        identifier set_left;
        identifier set_right;
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
    data_[0] = seqt_data{nil,0,0,0};
    atom_index_[0] = 0;
}

bool mem::advance_thread(thread& t, symbol c) {
    if(c == 0) return false;

    identifier atom_id = get_atom_id(c);

    seqt_data * tail_data = &data_[t.tail_.id_];
    seqt_data * following_data = &data_[t.following_.id_];

    std::set<identifier> visited;

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
        // it's repr's all the way down
        // for(identifier r = tail_data->repr; r != 0; r = data_[r].repr) {
        //     if(r == atom_id)
        //         return true;
        // }
        // let's keep it O(1) for now
        if(tail_data->repr == atom_id) {
            t.tail_ = tail_data->seq_next;
            return true;
        }
        return false;
        
    case set:
        // this time tail points to the set of found elements
        // assume it's a set of all atoms for now
        
        // look in the tail to see if we've already found it
        //TODO: we may need to check for circular loops

        for(;;) {
            if(visited.find(tail_data->id) != visited.end()) {
                // circular reference found!!
                throw std::logic_error("circular reference in set");
            }
            if(tail_data->type != set) {
                throw std::logic_error("set of non-set elements found!");
            }
                
            visited.insert(tail_data->id);

            if(tail_data->id < atom_id) {
                if(tail_data->set_left == 0)
                    break;
                tail_data = &data_[tail_data->set_left];
                continue;
            }
            else if(tail_data->id > atom_id) {
                if(tail_data->set_right == 0)
                    break;
                tail_data = &data_[tail_data->set_right];
                continue;
            }
            
            return false; // already found!
        }

        // we haven't found it, but tail_data points to where it should be
        // now look in following to see if it's even supposed to be part of this set

        for(;;) {
            if(visited.find(following_data->id) != visited.end()) {
                // circular reference found!!
                throw std::logic_error("circular reference in set");
            }
            if(tail_data->type != set) {
                throw std::logic_error("set of non-set elements found!");
            }
                
            visited.insert(following_data->id);

            if(following_data->id < atom_id) {
                if(following_data->set_left == 0)
                    break;
                following_data = &data_[following_data->set_left];
                continue;
            }
            else if(following_data->id > atom_id) {
                if(following_data->set_right == 0)
                    break;
                following_data = &data_[following_data->set_right];
                continue;
            }

            // we found it!
            break;
        }

        if(following_data->id == atom_id) {
            // insert into the tail
            if(tail_data->id < atom_id) {
                // insert right
            } else {
                // insert left
            }
            return true;
        }
        return false;
    }

    return false;
}



void mem::read_option1(symbol c) {
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