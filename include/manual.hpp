#ifndef __MANUAL_HPP__
#define __MANUAL_HPP__

#include <map>
#include <set>
#include <tuple>
#include <list>
#include <exception>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <algorithm>

using std::list;
using std::map;
using std::set;
using std::tuple;
using std::tie;
using std::get;
using std::vector;
using std::set_intersection;
using std::inserter;
using std::ostream;

struct seqt {
    enum node_type {
        atom = 0, sequence, collection
    };

    struct node;
    struct node_weight_less;

    typedef list<node*>::iterator sequence_iterator;
    typedef set<node*>::iterator collection_iterator;

    struct sequence_cursor {
        node * seq;
        size_t pos;
        sequence_iterator iter;

        sequence_cursor() 
            : seq(nullptr), pos(SIZE_MAX), iter() 
        { }
        sequence_cursor(node * n, size_t p) 
            : seq(n), pos(p)
        {
            if(seq == nullptr) return;
            
            iter = seq->seq.begin();
            for(size_t i = 0; i < p; i++)
                iter++;
        }
        sequence_cursor(node * n, size_t p, sequence_iterator i)
            : seq(n), pos(p), iter(i)
        { }
        sequence_cursor(sequence_cursor const & r) 
            : seq(r.seq), pos(r.pos), iter(r.iter)
        { }
        sequence_cursor(sequence_cursor && r) 
            : seq(r.seq), pos(r.pos), iter(r.iter)
        { }
        sequence_cursor & operator=(sequence_cursor const & r) {
            seq = r.seq;
            pos = r.pos;
            iter = r.iter;
            return *this;
        }
        sequence_cursor & operator=(sequence_cursor && r) {
            seq = r.seq; r.seq = nullptr;
            pos = r.pos; r.pos = 0;
            iter = r.iter; 
            return *this;
        }

        bool operator<(sequence_cursor r) const {
            if(seq < r.seq) return true;
            if(r.seq < seq) return false;
            return pos < r.pos;
        }

        static sequence_cursor lower(node * n) {
            return sequence_cursor{n, 0};
        }
    };
    
    // struct sequence_reference_set {
    //     map<node*, list<sequence_iterator>> refs;

    //     bool add(node* n, sequence_iterator i) {
    //         if(!refs.contains(n)) 
    //             refs[n] = list<sequence_iterator>();
            
    //         // slow but whatever
    //         for(auto j : refs[n]) {
    //             if(j == i) return false;
    //         }

    //         refs[n].push_back(i);
    //         return true;
    //     }

    //     bool remove(node * n, sequence_iterator i) {
    //         if(!refs.contains(n)) return false;

    //         list<sequence_iterator>::iterator pos = refs[n].begin();
    //         for(; pos != refs[n].end(); pos++)
    //             if(*pos == i) break;

    //         if(pos == refs[n].end())
    //             return false;

    //         refs[n].erase(pos);
    //     }

    //     bool contains(node * n, sequence_iterator j) {
    //         if(!refs.contains(n)) return false;

    //         for(auto i = refs[n].begin(); i != refs[n].end(); i++)
    //             if(*i == j) return true;
            
    //         return false;
    //     }

    //     struct const_iterator {
    //         map<node*, list<sequence_iterator>>::iterator map_position;
    //         map<node*, list<sequence_iterator>>::iterator map_end;
    //         list<sequence_iterator>::iterator list_position;

    //         const_iterator& operator++() {
    //             if(map_position == map_end) 
    //                 return *this;

    //             ++list_position;

    //             while(list_position == map_position->second.end()) {
    //                 ++map_position;
    //                 if(map_position == map_end)
    //                     return *this;

    //                 list_position = map_position->second.begin();
    //             }
    //         }
    //         const_iterator operator++(int) {
    //             const_iterator r = *this;
    //             this->operator++();
    //             return r;
    //         }
    //         tuple<node*, sequence_iterator> operator*() {
    //             if(map_position == map_end) 
    //                 return {nullptr, sequence_iterator()};
                
    //             return {map_position->first, *list_position};
    //         }
    //         bool operator==(const_iterator const & r) const {
    //             if(map_position != r.map_position) return false;
    //             if(map_position == map_end) return true;

    //             return list_position == r.list_position;
    //         }
    //         inline bool operator!=(const_iterator const & r) const {
    //             return !(*this == r);
    //         }
    //     };

    //     const_iterator begin() {
    //         auto b = refs.begin();
    //         if(b == refs.end())
    //             return const_iterator{b, refs.end()};
            
    //         return const_iterator{b, refs.end(), b->second.begin()};
    //     }
    //     const_iterator end() {
    //         const_iterator{refs.end(), refs.end()};
    //     }

        
    // };

    struct node {
        node_type type;

        size_t weight;

        uint32_t repr;
        list<node*> seq;
        set<node*> col;

        set<sequence_cursor> in_seq;
        map<node*, collection_iterator> in_col;

        node() 
            : type(atom), weight(1), repr(0) 
        { }
        node(uint32_t s) 
            : type(atom), weight(1), repr(s) 
        { }
        node(node_type t)
            : type(t), weight(1), repr(0)
        { }
        
        // don't call this unless you know what you are doing
        node(node const & n) 
            : type(n.type), weight(0), repr(n.repr), seq(n.seq), col(n.col), in_seq(), in_col()
        { /* now append something! */ }

        void add_sequence_reference(node * n, size_t pos, sequence_iterator i) {
            in_seq.insert({n, pos, i});
        }

        void remove_sequence_references(node * n) {
            auto i = in_seq.lower_bound(sequence_cursor::lower(n));
            if(i == in_seq.end()) 
                return;

            ++i;
            for(;i->seq == n; ++i)
                in_seq.erase(i);
        }

        bool expects(uint32_t s, list<node*> & chain, set<node*> & visited) {
            if(visited.contains(this)) 
                return false; // circular reference
            visited.insert(this);

            switch(type) {
            case atom:
                if(repr == s) {
                    chain.push_front(this);
                    return true;
                }
                return false;
            case sequence:
                // does the first node expect s?
                if(seq.front()->expects(s, chain, visited)) {
                    chain.push_front(this);
                    return true;
                }
                return false;
            case collection:
                // now we have to see if any of the nodes in the collection expect s
                // fallthrough
                break;
            }

            // we are a collection
            for(auto ci = col.begin(); ci != col.end(); ci++) {
                node * n = *ci;
                if(n->expects(s, chain, visited)) {
                    chain.push_front(this);
                    return true;
                }
            }

            return false;
        }

        void append(node * n) {
            switch(type) {
            case sequence:
                _append_sequence(n);
                break;
            case collection:
                _append_collection(n);
                break;
            case atom:
                // error out here
                break;
            }
            //TODO: add a logic_error or something
        }

        void _append_sequence(node * n) {
            // append n to the sequence
            size_t pos = seq.size();
            sequence_iterator i = seq.insert(seq.end(), n);
            // let n know where it is appended in our sequence
            n->add_sequence_reference(this, pos, i);
        }

        void _append_collection(node * n) {
            collection_iterator pos;
            bool success;

            tie(pos, success) = col.insert(n);
            if(!success) return;

            n->in_col.insert({this, pos});
        }
    };

    struct node_weight_less {
        bool operator()(node* l, node *r) const {
            if(r->weight < l->weight) return true;
            if(l->weight < r->weight) return false;
            return l < r; // just ensure they aren't the same node
        }
    };


    typedef set<node*,node_weight_less> node_set_type;
    typedef set<node*,node_weight_less>::iterator node_iterator;

    void remove_node(node *);
    node* make_atom(uint32_t s);
    node* make_seq(node*, node*);
    node* make_seq(sequence_iterator, sequence_iterator, node*);
    node* make_col(node*, node*);
    void read(uint32_t);
    void cleanup(size_t min_weight);
    bool is_duplicate(node_type, node * a, node * b);
    void read_node(
        node * this_one, list<node*> & todo,
        set<node*> & visited,
        set<sequence_cursor> & next_waiting_for,
        set<tuple<sequence_cursor, node*>> new_seq,
        set<tuple<node*, node*>> new_col);

    ostream & dump(ostream & os);

    ~seqt();

    node_set_type nodes;
    map<uint32_t,node*> atom_index;
    set<sequence_cursor> waiting_for;
    list<node*> completed;

};


seqt::node * seqt::make_atom(uint32_t s) {
    auto a = atom_index.find(s);
    if(a != atom_index.end())
        return a->second;

    node * n = new node(s);
    nodes.insert(n);
    atom_index[s] = n; // remember where this atom is

    return n;
}

seqt::node * seqt::make_seq(node * a, node * b) {
    node * r = new node(sequence);
    for(node * x : {a, b}) {
        if(x->type == sequence) {
            for(auto i = x->seq.begin(); i != x->seq.end(); i++) {
                r->append(*i);
            }
        } else {
            r->append(x);
        }
    }
    nodes.insert(r);
    return r;
}
seqt::node * seqt::make_seq(sequence_iterator b, sequence_iterator e, node * n) {
    node * r = new node(sequence);
    for(auto i = b; i != e; i++) {
        r->append(*i);
    }
    if(n->type == sequence) {
        for(auto i = n->seq.begin(); i != n->seq.end(); i++) {
            r->append(*i);
        }
    } else {
        r->append(n);
    }
    nodes.insert(r);
    return r;
}

seqt::node * seqt::make_col(node * a, node * b) {
    node * r = new node(collection);
    for(node * x : {a, b}) {
        if(x->type == collection) {
            for(auto i = x->col.begin(); i != x->col.end(); i++) {
                r->append(*i);
            }
        } else {
            r->append(x);
        }
    }
    nodes.insert(r);
    return r;
}

seqt::~seqt() {
    // gross but functional deletion
    auto temp = nodes; // copy it

    for(auto i = temp.begin(); i != temp.end(); i++) {
        remove_node(*i);
    }

    temp.clear();
}


void seqt::remove_node(node * r) {
    // unlink n from any sequences that reference it
    set<node*> to_remove;
    vector<node*> stack;
    stack.push_back(r);
    while(!stack.empty()) {
        node * cur = stack.back();
        stack.pop_back();

        if(to_remove.contains(cur))
            continue;
        to_remove.insert(cur);

        for(auto s : cur->in_seq) {
            stack.push_back(s.seq);
        }
    }

    for(node * n : to_remove) {
        // we checked for dependencies, can we just remove them?

        for(auto m : n->seq) {
            m->remove_sequence_references(n);
        }
        // unlink n from any collections that reference it
        for(auto j = n->in_col.begin(); j != n->in_col.end(); j++) {
            node * p = j->first;
            collection_iterator k = j->second;

            // should we replace it with a collection that n is in, or a collection of all the collections?
            p->col.erase(k);
        }
        n->in_col.clear();

        nodes.erase(n);
        delete n;
    }
}

void seqt::cleanup(size_t min_weight) {
    vector<node*> to_remove;

    for(node * n : nodes) {
        if(n->weight < min_weight) 
            to_remove.push_back(n);
    }

    for(node * n : to_remove) 
        remove_node(n);
}

void seqt::read(uint32_t s) {
    node * atom_node = make_atom(s);

    set<node*> visited;
    set<sequence_cursor> next_waiting_for;
    list<node*> todo;
    set<tuple<sequence_cursor, node*>> new_seq;
    set<tuple<node*, node*>> new_col;

    todo.push_back(atom_node);
    list<node*>::iterator cur = todo.begin();
    
    while(cur != todo.end()) {
        read_node(*cur, todo, visited, next_waiting_for, new_seq, new_col);
        cur++;
    }

    for(node * prev : completed) {
        for(node * c : todo) {
            if(c->type != atom) continue;

            switch(prev->type) {
            case atom:
            case collection:
                make_seq(prev, c);
                break;
            case sequence:
                make_seq(prev->seq.begin(), prev->seq.end(), c);
                break;
            }
            
        }
    }

    for(auto t : new_seq) {
        sequence_cursor cur;
        node * next;
        tie(cur, next) = t;
        // none of these should be duplicates ore we would have found them

        make_seq(cur.seq->seq.begin(), cur.iter, next);
    }

    for(auto t : new_col) {
        node * a, * b;
        tie(a, b) = t;
        // not sure how to check for duplicates, just add all for now
        make_col(a, b);
    }

    completed = todo;
    waiting_for = next_waiting_for;
}

void seqt::read_node(
    node * n, list<seqt::node*> & todo,
    set<seqt::node*> & visited, 
    set<sequence_cursor> & next_waiting_for,
    set<tuple<sequence_cursor, node*>> new_seq,
    set<tuple<node*, node*>> new_col) 
{
    if(visited.contains(n)) 
        return;
    visited.insert(n);

    n->weight++; // increment the weight of this node because we found it

    set<node*> potential_collections;
    set<sequence_cursor> potential_sequences;

    // get all the collections that this node is in
    for(auto p : n->in_col) {
        node * c;
        collection_iterator ci;
        tie(c, ci) = p;
        c->weight++; // increment this collection
        
        potential_collections.insert(c);
    }

    // now find any sequences that contain these collections, or the node itself
    // these reference the position of this node, so if it's at the beggining then we can ignore it
    for(auto s : n->in_seq)
        potential_sequences.insert(s);
        
    for(node * col : potential_collections) {   
        for(auto s : col->in_seq)
            potential_sequences.insert(s);
    }

    for(auto p : potential_sequences) {
        node * pn = p.seq;
        size_t pos = p.pos;
        sequence_iterator psi = p.iter, psn = psi;
        psn++;

        // this symbol is right at the beginning so add it to the next_waiting_for
        if(pn->seq.begin() == psi) {
            next_waiting_for.insert(sequence_cursor{pn, pos+1, psn});
            // next_waiting_for.add(pn, psn);
            continue;
        }

        // is this in our waiting for?
        if(waiting_for.contains(p)) {
            // add it to our list of matches
            if(psn == pn->seq.end())
                todo.push_back(pn);
            else
                next_waiting_for.insert(sequence_cursor{pn, pos+1, psn});

            waiting_for.erase(p);
            continue;
        } 
    }

    // whatever is left in waiting_for should be added as new sequences/cols
    for(auto sr : waiting_for) {
        new_seq.insert({sr, n});
        new_col.insert({sr.seq, n});
        for(auto c : potential_collections) {
            new_seq.insert({sr, c});
            new_col.insert({sr.seq, c});
        }
    }
}


// bool seqt::is_duplicate(seqt::node_type typ, seqt::node * a, seqt::node * b) {
//     if(typ == sequence) {
//         for(auto st : a->in_seq) {

//             // is si the at the begining?
//             if(si != s->seq.begin())
//                 continue;
        
//             // does b come after a?
//             si++;
//             if(si == s->seq.end() || *si != b) 
//                 continue;

//             // is that the end?
//             si++;
//             if(si != s->seq.end())
//                 continue;

//             // if we are still here then we already have a sequence of a and b
//             return true;
//         }
//     } else {
//         // collection
//         for(auto co : a->in_col) {
//             node * c;
//             collection_iterator ci;
//             tie(c, ci) = co;

//             // does this collection only have two things?
//             if(c->col.size() != 2) 
//                 continue;

//             // are the two things a and b?
//             if(!c->col.contains(b))
//                 continue;

//             return true;
//         }
//     }

//     return false;
// }


ostream & seqt::dump(ostream & os) {
    uint32_t id = 0;

    map<node*, uint32_t> ids;

    auto write_char = [&os](uint32_t c) {
        if(c >= 33 && c <= 126) {
            os << (char)c;
            return;
        }
        os << "\\u" << std::hex << std::setfill('0') << std::setw(4) << (uint16_t)c << std::dec;
    };

    os << "[\n";
    for(node * n : nodes) {
        ids[n] = id++;
        bool first = true;

        os << "{\n";
        os << "\t\"id\": " << ids[n] << "\n";
        os << "\t\"weight\": " << n->weight << "\n";
        switch(n->type) {
        case atom:
            os << "\t\"atom\": \"";
            write_char(n->repr);
            os << "\"\n";
            break;
        case sequence:
            os << "\t\"seq\": [";
            first = true;
            for(node * s : n->seq) {
                if(first) first = false;
                else os << ",";
                if(!ids.contains(s))
                    ids[s] = id++;
                
                os << ids[s];
            }
            os << "]\n";
            break;
        case collection:
            os << "\t\"col\": [";
            first = true;
            for(node * c : n->col) {
                if(first) first = false;
                else os << ",";

                if(!ids.contains(c))
                    ids[c] = id++;
                
                os << ids[c];
            }
            os << "]\n";
            break;
        }
        os << "}\n";
    }
    os << "]\n";
    return os;
}



#endif