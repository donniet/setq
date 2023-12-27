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

using std::list;
using std::map;
using std::set;
using std::tuple;
using std::tie;
using std::get;
using std::vector;
using std::set_intersection;
using std::inserter;

struct seqt {
    enum node_type {
        atom = 0, sequence, collection
    };

    struct node;
    struct sequence_reference;
    struct node_weight_less;

    typedef vector<node*>::iterator sequence_iterator;
    typedef set<node*>::iterator collection_iterator;

    struct sequence_reference
        : public tuple<node*, sequence_iterator>
    {
        bool operator<(sequence_reference const & r) const {
            return get<0>(*this) < get<0>(r);
                return true;
            if(get<0>(r) > get<0>(*this))
                return false;
            
            // if the two nodes are equal, we will just compare the iterators
            sequence_iterator i = get<1>(*this);
            sequence_iterator j = get<1>(r);

            return i < j; // we can do this because these are vector iterators
        }

        sequence_reference(node* n, sequence_iterator i) 
            : tuple<node*, sequence_iterator>({n, i}) 
        { }
    };

    struct node {
        node_type type;

        size_t weight;

        uint32_t repr;
        vector<node*> seq;
        set<node*> col;

        set<sequence_reference> in_seq;
        map<node*, collection_iterator> in_col;

        node() 
            : type(atom), weight(0), repr(0) 
        { }
        node(uint32_t s) 
            : type(atom), weight(0), repr(s) 
        { }
        
        // don't call this unless you know what you are doing
        node(node const & n) 
            : type(n.type), weight(0), repr(n.repr), seq(n.seq), col(n.col), in_seq(), in_col()
        { /* now append something! */ }

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
            sequence_iterator pos = seq.insert(seq.end(), n);
            // let n know where it is appended in our sequence
            n->in_seq.insert(sequence_reference(this, pos));
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

    void remove_node(node_iterator i);
    node_iterator make_atom(uint32_t s);
    node_iterator make_seq(node_iterator, node_iterator);
    node_iterator make_col(node_iterator, node_iterator);
    void read(uint32_t);
    bool is_duplicate(node_type, node * a, node * b);
    void read_node(
        vector<node*>::iterator this_one, vector<node*> & todo, 
        set<node*> & visited, set<sequence_reference> & next_waiting_for,
        set<tuple<node_type, node*, node*>> & new_nodes);

    ~seqt();

    node_set_type nodes;
    map<uint32_t,node_iterator> atom_index;
    set<sequence_reference> waiting_for;

};

seqt::node_iterator seqt::make_atom(uint32_t s) {
    auto a = atom_index.find(s);
    if(a != atom_index.end())
        return a->second;

    node * n = new node(s);
    bool _; // insertion should always succeed
    node_iterator i;
    tie(i, _) = nodes.insert(n);
    atom_index[s] = i; // remember where this atom is

    return i;
}

seqt::~seqt() {
    // gross but functional deletion
    for(auto i = nodes.begin(); i != nodes.end(); i++) {
        remove_node(i);
    }
}


void seqt::remove_node(node_iterator i) {
    node * n = *i;

    // unlink n from any sequences that reference it
    for(auto j = n->in_seq.begin(); j != n->in_seq.end(); j++) {
        node * p;
        sequence_iterator k;
        tie(p, k) = *j;

        // should we replace it with a collection of all the collections that n is in maybe?
        // so many possibilities...
        p->seq.erase(k);
    }

    // unlink n from any collections that reference it
    for(auto j = n->in_col.begin(); j != n->in_col.end(); j++) {
        node * p = j->first;
        collection_iterator k = j->second;

        // should we replace it with a collection that n is in, or a collection of all the collections?
        p->col.erase(k);
    }

    // now remove n from our list
    nodes.erase(i);
    delete n;
}


void seqt::read(uint32_t s) {
    node_iterator ni = make_atom(s);
    node * atom_node = *ni;

    set<node*> visited;
    set<sequence_reference> next_waiting_for;
    vector<node*> todo;
    set<tuple<node_type, node*, node*>> new_nodes;

    todo.push_back(atom_node);
    vector<node*>::iterator cur = todo.begin();
    
    while(cur != todo.end()) {
        read_node(todo.begin(), todo, visited, next_waiting_for, new_nodes);
        cur++;
    }

    // go through all the new potential nodes and try and remove any duplicates before creating them
    for(auto t : new_nodes) {
        node_type typ;
        node * a, * b;
        tie(typ, a, b) = t;

        if(!is_duplicate(typ, a, b)) {
            node * c = new node();
            c->type = typ;
            c->append(a);
            c->append(b);

            nodes.insert(c);  // track our new node
        }
    }

    waiting_for = next_waiting_for;
}

bool seqt::is_duplicate(seqt::node_type typ, seqt::node * a, seqt::node * b) {
    if(typ == sequence) {
        for(auto st : a->in_seq) {
            node * s;
            sequence_iterator si;
            tie(s, si) = st;

            // is si the at the begining?
            if(si != s->seq.begin())
                continue;
        
            // does b come after a?
            si++;
            if(si == s->seq.end() || *si != b) 
                continue;

            // is that the end?
            si++;
            if(si != s->seq.end())
                continue;

            // if we are still here then we already have a sequence of a and b
            return true;
        }
    } else {
        // collection
        for(auto co : a->in_col) {
            node * c;
            collection_iterator ci;
            tie(c, ci) = co;

            // does this collection only have two things?
            if(c->col.size() != 2) 
                continue;

            // are the two things a and b?
            if(!c->col.contains(b))
                continue;

            return true;
        }
    }

    return false;
}

void seqt::read_node(
    vector<seqt::node*>::iterator this_one, vector<seqt::node*> & todo, 
    set<seqt::node*> & visited, set<seqt::sequence_reference> & next_waiting_for,
    set<tuple<seqt::node_type, seqt::node*, seqt::node*>> & new_nodes) 
{
    node * n = *this_one;

    if(visited.contains(n)) 
        return;
    visited.insert(n);

    n->weight++; // increment the weight of this node because we found it


    set<node*> potential_collections;
    set<sequence_reference> potential_sequences;

    // get all the collections that this node is in
    for(auto p : n->in_col) {
        node * c;
        collection_iterator ci;
        tie(c, ci) = p;
        c->weight++; // increment this collection
        
        potential_collections.insert(c);
    }

    auto find_sequences = [&](node * n) -> bool {
        bool inserted = false;
        for(auto seqj = n->in_seq.begin(); seqj != n->in_seq.end(); seqj++) {
            node * s;
            sequence_iterator si;
            tie(s, si) = *seqj;

            set<sequence_reference>::iterator _;
            bool ins;
            
            tie(_, ins) = potential_sequences.insert(sequence_reference(s, si));
            inserted = (inserted || ins);
        }
        return inserted;
    };

    // now find any sequences that contain these collections, or the node itself
    bool found = find_sequences(n);
    for(node * col : potential_collections) {
        bool found_col = find_sequences(col);
        found = (found || found_col);
    }

    // are any of these sequences ones we are waiting for?
    set<sequence_reference> matches;

    for(auto p : potential_sequences) {
        node * pn;
        sequence_iterator psi;
        tie(pn, psi) = p;

        // is this in our waiting for?
        if(waiting_for.contains(p)) {
            // add it to our list of matches
            matches.insert(p);
            continue;
        }
        // otherwise, maybe we should be tracking these?
        for(auto w : waiting_for) {
            node * wn;
            sequence_iterator wsi;
            tie(wn, wsi) = w;

            // perhaps we are missing a collection made from the two expected nodes
            new_nodes.insert({collection, *psi, *wsi});
            // or perhaps we are missing a sequence with this node and the node we are waiting for
            new_nodes.insert({sequence, n, *wsi});
            // TODO: think of other things?
        }
    }

    set<sequence_reference> matches_next;
    // now move all the matched ones forward one and check for end
    for(auto p : matches) {
        node * n;
        sequence_iterator si;
        tie(n, si) = p;

        ++si;
        if(si == n->seq.end()) {
            // we've completed this sequence and now need to look for 
            // add this node to or to_do list
            todo.push_back(n);
        } else {
            // otherwise add this to the potential sequences we are still waiting on
            next_waiting_for.insert(sequence_reference(n, si));
        }
    }
}



#endif