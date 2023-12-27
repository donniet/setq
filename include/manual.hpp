#ifndef __MANUAL_HPP__
#define __MANUAL_HPP__

#include <map>
#include <set>
#include <tuple>
#include <list>

using std::list;
using std::map;
using std::set;
using std::tuple;
using std::tie;
using std::get;

struct seqt {
    typedef enum {
        atom = 0, sequence, collection
    } node_type;

    struct node;

    std::list<node*> nodes;

    typedef std::list<node*>::iterator node_iterator;

    void remove_node(node_iterator i);
};

struct seqt::node {
    node_type type;

    size_t weight;

    uint32_t repr;
    list<node*> seq;
    set<node*> col;

    typedef list<node*>::iterator sequence_iterator;
    typedef set<node*>::iterator collection_iterator;

    set<tuple<node*, sequence_iterator>> in_seq;
    map<node*, collection_iterator> in_col;

    node() : type(atom), weight(0), repr(0) { }

    void append(node * n) {
        switch(type) {
        case sequence:
            _append_sequence(n);
            break;
        case collection:
            _append_collection(n);
            break;
        }
        //TODO: add a logic_error or something
    }

    void _append_sequence(node * n) {
        // append n to the sequence
        sequence_iterator pos = seq.insert(seq.end(), n);
        // let n know where it is appended in our sequence
        n->in_seq.insert({this, pos});
    }

    void _append_collection(node * n) {
        collection_iterator pos;
        bool success;

        tie(pos, success) = col.insert(n);
        if(!success) return;

        n->in_col.insert({this, pos});
    }
};

void seqt::remove_node(node_iterator i) {
    node * n = *i;

    // unlink n from any sequences that reference it
    for(auto j = n->in_seq.begin(); j != n->in_seq.end(); j++) {
        node * p;
        node::sequence_iterator k;
        tie(p, k) = *j;

        // should we replace it with a collection of all the collections that n is in maybe?
        // so many possibilities...
        p->seq.erase(k);
    }

    // unlink n from any collections that reference it
    for(auto j = n->in_col.begin(); j != n->in_col.end(); j++) {
        node * p = j->first;
        node::collection_iterator k = j->second;

        // should we replace it with a collection that n is in, or a collection of all the collections?
        p->col.erase(k);
    }

    // now remove n from our list
    nodes.erase(i);
}

#endif