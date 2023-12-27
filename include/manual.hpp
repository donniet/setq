#ifndef __MANUAL_HPP__
#define __MANUAL_HPP__

#include <map>
#include <set>
#include <list>

using std::list;
using std::map;
using std::set;

struct seqt {
    typedef uint32_t ident;

    typedef enum {
        atom = 0, sequence, collection
    } node_type;

    struct node;
};

struct seqt::node {
    ident id;
    node_type type;
    
    uint32_t repr;
    list<ident> seq;
    set<ident> col;

    map<ident, list<ident>::const_iterator> in_seq;
    map<ident, list<ident>::const_iterator> in_col;
};

#endif