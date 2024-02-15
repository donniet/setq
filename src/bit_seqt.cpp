#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <execution>
#include <functional>
#include <numeric>
#include <set>
#include <cmath>
#include <string>
#include <sstream>
#include <memory>

// #ifndef DEBUG
#  if ( __clang__ && __APPLE__ )
// disable parallel exeuction for now
#  else
#    define EXECUTION_PARALLEL 1
#  endif
// #endif

#if EXECUTION_PARALLEL == 1
#  define PAR_UNSEQ std::execution::par_unseq,
#else
#  define PAR_UNSEQ
#endif


using namespace std;

struct bit_seqt {
    struct node;
    struct node_ptr;
    struct node_pair;

    static constexpr float statistical_significance = 9.; // standard deviations
    static constexpr size_t max_allocation = 16e9; // 16GB

    size_t capacity;
    vector<node> nodes;
    vector<node_pair> node_index;
    vector<node_pair> suggested;
    vector<size_t> active_scrap;
    vector<size_t> preactive_scrap;


    vector<size_t> nodes_active;
    vector<size_t> nodes_active_scan;
    vector<size_t> activation;
    vector<size_t> activation_scan;
    vector<float> charge_buffer;


    static float pvalue(size_t first, size_t second, size_t count_sequence);

    node_ptr ptr(size_t offset);
    node_ptr ptr(node const & n);

    void grow_capacity(size_t new_capacity = 0);

    void save(ostream & os);
    void read(bool bit);
    void prune(size_t max_nodes);
    bool pair_exists(node_ptr first, node_ptr second);
    node_ptr ancestor(node_ptr first, node_ptr second);
    node create_relationship(node_ptr first, node_ptr second);
    bool is_novel_suggestion(node_ptr first, node_ptr second);
    void pack_suggestions(node & n, node_ptr next, vector<node>::iterator output_begin);
    void pack_suggestion_index(node & n, node_ptr next, vector<node_pair>::iterator output_begin);

    bit_seqt(size_t initial_size = 1024);
};

// number of standard deviations from the expected value ab is
float sequence_stddevs(size_t a, size_t b, size_t ab) {
    if(ab > a || ab > b || a + b == 0)
        return 0;

    float total = (float)(a + b);
    float pa = (float)a / total;
    float pb = (float)b / total;
    float pab = pa * pb;
    float expected_value = total * pab;
    float variation = total * pab * (1-pab);
    float stddev = sqrt(variation);

    return ((float)ab - expected_value) / stddev;
}

float sequence_pvalue(size_t a, size_t b, size_t ab) {
    if(a == 0 && b == 0) {
        if(ab > 0) return 0.;
        return 0.5;
    }
    float total = (float)(a + b);
    float pa = (float)a / total;
    float pb = (float)b / total;
    float pab = pa * pb;
    float expected_value = total * pab;
    float variation = total * pab * (1-pab);

    float x = ((float)ab - expected_value) * sqrt(2 / variation);

    return 0.5 * erfc(x);
}

float bit_seqt::pvalue(size_t first, size_t second, size_t sequence) {
    return sequence_pvalue(first, second, sequence);
}


struct bit_seqt::node_ptr {
private:
    vector<node> * nodes_;
    size_t offset_;
public:
    node_ptr() : nodes_(nullptr), offset_(0) { }
    static node_ptr null() { return node_ptr(); }

    node_ptr(vector<node> & n, size_t o) : nodes_(&n), offset_(o) { }
    node_ptr(node_ptr const & r) : nodes_(r.nodes_), offset_(r.offset_) { }
    node_ptr(node_ptr && r) : nodes_(r.nodes_), offset_(r.offset_) {
        r.nodes_ = nullptr;
        r.offset_ = 0;
    }

    node_ptr & operator=(node_ptr const & r) {
        nodes_ = r.nodes_;
        offset_ = r.offset_;
        return *this;
    }
    node_ptr & operator=(node_ptr && r) {
        nodes_ = r.nodes_;
        offset_ = r.offset_;
        r.nodes_ = nullptr;
        r.offset_ = 0;
        return *this;
    }
    
    bool operator==(node_ptr const & r) const {
        if(nodes_ == nullptr && r.nodes_ == nullptr) return true;
        if(nodes_ == nullptr || r.nodes_ == nullptr) return false;
        return nodes_ == r.nodes_ && offset_ == r.offset_;
    }
    bool operator!=(node_ptr const & r) const {
        return !(*this == r);
    }
    bool operator==(void const * r) const {
        return (node*)(*this) == r;
    }
    bool operator!=(void const * r) const {
        return (node*)(*this) != r;
    }

    void reset() {
        nodes_ = nullptr;
        offset_ = 0;
    }

    bool is_null() const;
    bool operator<(node_ptr const & r) const;
    bool operator<=(node_ptr const & r) const;

    operator node*() const;
    operator size_t() const {
        return offset_;
    }
    size_t offset() const {
        return offset_;
    }

    node & operator*() const;
    node * operator->() const;
    node_ptr & operator++();
    node_ptr & operator--();
    node_ptr operator++(int);
    node_ptr operator--(int);
};


struct bit_seqt::node {
    enum node_operation {
        atom, sequence, relationship
    } op;
    size_t count;
    float charge;

    node_ptr first;
    node_ptr second;
    node_ptr parent;

    node() : op(atom), count(0), charge(0), first(), second(), parent() { } 
    node(node_operation o, node_ptr f, node_ptr s) :
        op(o), count(0), charge(0), first(f), second(s), parent()
    { }
    node(node const & r) : 
        op(r.op), count(r.count), charge(r.charge), first(r.first), second(r.second), parent(r.parent)
    { }
    node(node && r) : 
        op(r.op), count(r.count), charge(r.charge), first(r.first), second(r.second), parent(r.parent)
    { 
        r.op = atom;
        r.count = 0;
        r.charge = 0;
        r.first.reset();
        r.second.reset();
        r.parent.reset();
    }

    node & operator=(node const & r) {
        op = r.op;
        count = r.count;
        charge = r.charge;
        first = r.first;
        second = r.second;
        parent = r.parent;
        return *this;
    }
    node & operator=(node && r) {
        op = r.op;          r.op = atom;
        count = r.count;    r.count = 0;
        charge = r.charge;  r.charge = 0;
        first = r.first;    r.first.reset();
        second = r.second;  r.second.reset();
        parent = r.parent;  r.parent.reset();
        return *this;
    }

    void save(bit_seqt * s, ostream & os) const;
    void load(bit_seqt * s, istream & is);

    float significance() const;
    float activate();
    bool is_active() const;
};

bool bit_seqt::node_ptr::operator<(node_ptr const & r) const {
    // return &nodes_->operator[](0) + offset_ < &r.nodes_->operator[](0) + r.offset_;
    // it's fine to just use the nodes_ pointer as the reference for a less than comparison
    // and it avoids a null check
    return nodes_ + offset_ < r.nodes_ + r.offset_;
    
}
bool bit_seqt::node_ptr::operator<=(node_ptr const & r) const {
    // return &nodes_->operator[](0) + offset_ <= &r.nodes_->operator[](0) + r.offset_;
    // it's fine to just use the nodes_ pointer as the reference for a less than comparison
    // and it avoids a null check
    return nodes_ + offset_ <= r.nodes_ + r.offset_;
}
bit_seqt::node & bit_seqt::node_ptr::operator*() const {
    return nodes_->operator[](offset_);
}
bit_seqt::node * bit_seqt::node_ptr::operator->() const {
    return &nodes_->operator[](offset_);
}
bit_seqt::node_ptr::operator node*() const {
    // we have to do a null check here I believe
    return is_null() ? nullptr : &nodes_->operator[](0) + offset_;
}
bool bit_seqt::node_ptr::is_null() const {
    return nodes_ == nullptr || offset_ >= nodes_->size();
}
bit_seqt::node_ptr & bit_seqt::node_ptr::operator++() {
    if(offset_ < nodes_->size()) 
        offset_++;

    return *this;
}
bit_seqt::node_ptr & bit_seqt::node_ptr::operator--() {
    if(offset_ > 0)
        offset_--;

    return *this;
}
bit_seqt::node_ptr bit_seqt::node_ptr::operator++(int) {
    node_ptr r = *this;
    operator++();
    return r;
}
bit_seqt::node_ptr bit_seqt::node_ptr::operator--(int) {
    node_ptr r = *this;
    operator--();
    return r;
}

float bit_seqt::node::significance() const {

    if(count < 10) return 0;

    switch(op) {
    case atom:
        return numeric_limits<float>::max(); // atoms are infinitely significant
    case sequence:
        return sequence_stddevs(first->count, second->count, count);
    case relationship:
        return sequence_stddevs(first->first->count, first->second->count, first->count + second->count);
    }
    return 0.;
}

bool bit_seqt::node::is_active() const {
    switch(op) {
    case atom:
        // atoms are activated by the read function
        if(charge >= 1.) 
            return true;

        break;
    case sequence:
        if((first == second && first->charge >= 1.5) ||
           (second->charge >= 1. && first->charge >= 0.5 && first->charge < 1.)) 
            return true;

        break;
    case relationship:
        if(first->first->charge >= 1. || first->second->charge >= 1.)
            return true;
        
        break;
    }

    return false;
}

float bit_seqt::node::activate() {
    // take the charge from children
    // TODO: should this be all the charge
    switch(op) {
    case atom:
        // atoms are activated by the read function
        break;
    case sequence:
        if((first == second && first->charge >= 1.5) ||
           (second->charge >= 1. && first->charge >= 0.5 && first->charge < 1.)) 
        {
            count++;
            // TODO: I don't know if we should divide by 2 here
            return charge + 1.;
        }
        break;
    case relationship:
        if(first->first->charge >= 1. || first->second->charge >= 1.) {
            count++;
            return charge + 1.;
        }
        break;
    }

    return charge;
}



struct bit_seqt::node_pair {
    node_ptr first;
    node_ptr second;
    node_ptr parent;

    node_pair() : first(), second(), parent() {}
    node_pair(node_ptr f, node_ptr s, node_ptr p) : parent(p) {
        if(f < s) {
            first = f; 
            second = s;
        } else {
            first = s;
            second = f;
        }
    }

    bool operator<(node_pair const & r) const {
        if(first < r.first) return true;
        if(r.first < first) return false;
        return second < r.second;
    }
};


bit_seqt::node_ptr bit_seqt::ptr(size_t offset) {
    return node_ptr(nodes, offset);
}
bit_seqt::node_ptr bit_seqt::ptr(node const & n) {
    return node_ptr(nodes, &n - &nodes[0]);
}

void bit_seqt::prune(size_t max_nodes) {
    // // first calculate a usefulness of each node
    // vector<float> use(nodes.size());
    // scrap_nodes.resize(nodes.size());

    // //   get an iota equal to the number of nodes
    // for_each(PAR_UNSEQ nodes.begin(), nodes.end(), [&](node const & n) {
    //     size_t dex = &n - &nodes[0];

    //     // count the dependencies;
    //     scrap_nodes[dex] = dex;
    // });

    // //   save the usefulness of each node
    // for_each(PAR_UNSEQ nodes.begin(), nodes.end(), [&](node const & n) {
    //     size_t dex = &n - &nodes[0];

    //     switch(n.op) {
    //     case node::atom:
    //         use[dex] = numeric_limits<float>::max(); // infinitely useful
    //         break;
    //     case node::sequence:
    //         use[dex] = sequence_stddevs(n.first->count, n.second->count, n.count);
    //         break;
    //     case node::relationship:
    //         use[dex] = -sequence_stddevs(n.first->first->count, n.first->second->count, 
    //                                      n.first->count + n.second->count);
    //         break;
    //     }
    // });

    // flow the usefulness values down through child nodes
    // (if a parent node is useful, then the child's usefulness should be greater)

    // sort the usefulness values and find the cutoff

    // remove the least useful nodes

    // could we combine or restructure the nodes also?  so far this code creates essentially
    // a linked list of bits.  
}

void bit_seqt::read(bool bit) {
    using namespace std::placeholders;

    size_t current_size, index_size, active_count, preactive_count;

    
    // charge_buffer.resize(nodes.size());

    // add one to the node representing the bit read
    node_ptr cur = ptr(bit ? 1 : 0);
    cur->charge += 1.; // charge the atom
    cur->count++;

    nodes_active.resize(nodes.size());
    nodes_active_scan.resize(nodes.size());
    charge_buffer.resize(nodes.size());


    // this will keep track of the nodes that changed
    fill(PAR_UNSEQ nodes_active.begin(), nodes_active.end(), 0);
    nodes_active[cur.offset()] = 1;
    active_count = 1;

    // copy the charge into the charge buffer
    transform(PAR_UNSEQ nodes.begin(), nodes.end(), charge_buffer.begin(),
        [&](node const & n) { return n.charge; });

    // activate until all nodes to be activated have been
    for(;;) {
        // activate charge2 using the charge variables
        // TODO: filter out the nodes that already changed
        for_each(PAR_UNSEQ nodes.begin(), nodes.end(),
        [&](node & n) {
            size_t dex = &n - &nodes[0];

            // cout << "nodes[" << dex << "].charge = " << n.charge << endl;

            // don't re-activate an already activated node
            if(nodes_active[dex])
                return;

            if(n.is_active()) {
                n.count++;
                charge_buffer[dex] += 1.;
                nodes_active[dex] = 1;
            }
        });

        // unbuffer the new charge2 values back into charge
        // TODO: should the charge of child nodes go negative or somehow
        //   signal that a parent has picked up this node?
        for_each(PAR_UNSEQ charge_buffer.begin(), charge_buffer.end(), 
        [&](float const & charge) {
            size_t dex = &charge - &charge_buffer[0];
            if(nodes_active[dex])
                nodes[dex].charge = charge;
        });

        // count all the changed nodes
        inclusive_scan(PAR_UNSEQ nodes_active.begin(), nodes_active.end(), 
            nodes_active_scan.begin());

        // if none have changed we are done
        if(nodes_active_scan.back() == active_count)
            break;

        // othewise save the new count that has changed, and keep going
        active_count = nodes_active_scan.back();
    }

    // rotate active index;
    active_scrap.swap(preactive_scrap); 

    // ensure enough space in the scrap_index for the active node indices.
    active_scrap.resize(active_count);

    // pack the active node indexes into node_index
    for_each(PAR_UNSEQ nodes_active_scan.begin(), nodes_active_scan.end(),
    [&](size_t const & i) {
        size_t dex = &i - &nodes_active_scan[0];
        
        if(i == 0 || (dex > 0 && i == nodes_active_scan[dex-1]))
            return;

        active_scrap[i-1] = dex;
    });

    // minus 1 because this is an inclusive scan
    size_t new_node_index_begin = node_index.size();
    size_t new_node_begin = nodes.size();

    // look for sequences of the same node (we'll use the nodes_active vector as scrap)
    activation.resize(active_count);
    activation_scan.resize(active_count);

    for_each(PAR_UNSEQ active_scrap.begin(), active_scrap.end(), 
    [&](size_t const & active_i) {
        size_t dex = &active_i - &active_scrap[0];

        // cout << "nodes[" << dex << "].charge = " << nodes[dex].charge << endl;

        // is this node in sequence with itself?
        // and is that node
        if( nodes[dex].charge >= 1.5 && 
            // nodes[dex].count > 1000 &&
           !pair_exists(ptr(dex), ptr(dex)) &&
            abs(nodes[dex].significance()) > statistical_significance) 
        {
            activation[dex] = 1;
            return;
        }
        
        activation[dex] = 0;
    });

    // how many self-sequences will we create?
    inclusive_scan(PAR_UNSEQ activation.begin(), activation.end(),
        activation_scan.begin());
    active_count = activation_scan.back();

    // pack the new nodes and add their indexes
    if(active_count > 0) {
        if(2*(nodes.size() + active_count) >= capacity) {
            grow_capacity(4*(nodes.size() + active_count));
        }

        node_index.resize(node_index.size() + active_count);
        nodes.resize(nodes.size() + active_count);

        for_each(PAR_UNSEQ activation_scan.begin(), activation_scan.end(), 
        [&](size_t const & i) {
            size_t dex = &i - &activation_scan[0];

            if(i == 0 || (dex > 0 && i == activation_scan[dex-1]))
                return;

            node_ptr c = ptr(dex);
            
            nodes[new_node_begin + (i - 1)] = node(node::sequence, c, c);
            node_index[new_node_index_begin + (i - 1)] = node_pair(c, c, ptr(new_node_begin + (i - 1)));
        });

        new_node_index_begin = node_index.size();
        new_node_begin = nodes.size();

        sort(PAR_UNSEQ node_index.begin(), node_index.end());
    }


    // now we will create any new sequences that do not exist
    // but I don't know how to do this without it being NxM, so 
    // just reserviing NxM for now and we'll optimize later...
    // TODO: Maybe this could be probabilistic?  choose N potential pairs
    //  randomly and check those.  N can be chosen based on a certain bitrate
    //  we are aiming at.
    activation.resize(active_scrap.size() * preactive_scrap.size());
    activation_scan.resize(active_scrap.size() * preactive_scrap.size());

    if(activation.size() > 0) {
        for_each(PAR_UNSEQ activation.begin(), activation.end(),
        [&](size_t & c) {
            size_t i = &c - &activation[0];
            size_t active_j = i / preactive_scrap.size();
            size_t preact_k = i % preactive_scrap.size();

            c = 0;

            if(active_scrap[active_j] <= preactive_scrap[preact_k])
                return;

            node_ptr a(nodes, active_scrap[active_j]), 
                     p(nodes, preactive_scrap[preact_k]);

            // if(a->count < 1000 && p->count < 1000)
            //     return;
            
            if(abs(p->significance()) < statistical_significance)
                return;

            if(abs(a->significance()) < statistical_significance)
                return;
            
            // TODO: this function is not a great parallel function...
            node_ptr g = ancestor(a, p);
            if(!g.is_null())
                return;

            // now look to see if a or p contain the other
            // find a in p:

            c = 1;
        });

        inclusive_scan(PAR_UNSEQ activation.begin(), activation.end(), 
            activation_scan.begin());
        active_count = activation_scan.back();

        if(2*(nodes.size() + 3*active_count) >= capacity) {
            grow_capacity(4*(nodes.size() + 3*active_count));
        }

        node_index.resize(node_index.size() + active_count);
        nodes.resize(nodes.size() + 3 * active_count);

        for_each(PAR_UNSEQ activation_scan.begin(), activation_scan.end(),
        [&](size_t const & c) {
            size_t i = &c - &activation_scan[0];
            size_t active_j = i / preactive_scrap.size();
            size_t preact_k = i % preactive_scrap.size();

            if(c == 0 || (i > 0 && c == activation_scan[i-1]))
                return;

            node_ptr a(nodes, active_scrap[active_j]), 
                     p(nodes, preactive_scrap[preact_k]);

            // TODO: find the right insertion point to make relationships a valid binary tree

            // calculate the offset to this new node (each new node actually creates 3 node entries)
            size_t off = 3*(c-1);
            
            node & r = nodes[new_node_begin + off] = node(node::sequence, p, a);
            r.count = 1;
            r.charge = 0;

            node & s = nodes[new_node_begin + off + 1] = node(node::sequence, a, p);
            s.count = 0;
            s.charge = 0;

            // now we have to find the right insertion point for this new node:
            node & t = nodes[new_node_begin + off + 2] = create_relationship(a, p);

            // node & t = nodes[new_node_begin + off + 2] = node(node::relationship, ptr(r), ptr(s));
            // t.count = a->count + p->count;
            // t.charge = 0;

            // now add this node to the node_index
            off = (c-1);
            node_index[new_node_index_begin + off] = node_pair(p, a, ptr(t));
        });

        new_node_index_begin = node_index.size();
        new_node_begin = nodes.size();

        sort(PAR_UNSEQ node_index.begin(), node_index.end());
    }

    // now divide all the charges by 2 to get ready for the next run
    for_each(PAR_UNSEQ nodes.begin(), nodes.end(),
    [&](node & n) {
        n.charge /= 2.;
    });

}


void bit_seqt::save(ostream & os) {
    for(node const & n : nodes) {
        n.save(this, os);
    }
}

bool bit_seqt::pair_exists(node_ptr first, node_ptr second) {
    return binary_search(node_index.begin(), node_index.end(), node_pair(first, second, node_ptr::null()));
}

bit_seqt::node bit_seqt::create_relationship(node_ptr a, node_ptr b) {
    for(;;) {
        if(abs(a->significance()) > abs(b->significance()))
            swap(a, b);
        
        if(!b->parent.is_null())
            b = b->parent;

        if(!a->parent.is_null())
            a = a->parent;

        // both parents are null
        break;
    }


    node ret;
    ret.op = node::relationship;
    ret.first = a;
    ret.second = b;
    ret.count = a->count + b->count;
    ret.charge = 0.;
    return ret;
}

bit_seqt::node_ptr bit_seqt::ancestor(node_ptr first, node_ptr second) {
    vector<node_ptr> first_parents;
    for(node_ptr p = first; !p.is_null(); p = p->parent) {
        first_parents.push_back(p);
    }

    sort(PAR_UNSEQ first_parents.begin(), first_parents.end());

    for(node_ptr q = second; !q.is_null(); q = q->parent) {
        if(binary_search(first_parents.begin(), first_parents.end(), q))
            return q;
    }

    return node_ptr::null();
}

void bit_seqt::grow_capacity(size_t new_capacity) {
    if(new_capacity == 0)
        new_capacity = 2 * capacity;

    if((sizeof(nodes.front()) + sizeof(node_index.front()) + sizeof(suggested.front()) + sizeof(active_scrap.front()) + 
        sizeof(preactive_scrap.front()) + sizeof(nodes_active.front()) + sizeof(nodes_active_scan.front()) + sizeof(charge_buffer.front())) 
            * new_capacity + 
       (sizeof(activation.front()) + sizeof(activation_scan.front())) 
            * (new_capacity * new_capacity) 
        > max_allocation) 
    {
        throw logic_error("bad_alloc");
    }

    capacity = new_capacity;

    nodes.reserve(capacity);
    
    node_index.reserve(capacity);
    suggested.reserve(capacity);
    active_scrap.reserve(capacity);
    preactive_scrap.reserve(capacity);
    activation.reserve(capacity * capacity);
    activation_scan.reserve(capacity * capacity);

    nodes_active.reserve(capacity);
    nodes_active_scan.reserve(capacity);
    charge_buffer.reserve(capacity);
}

bit_seqt::bit_seqt(size_t initial_size) 
    : capacity(initial_size)
{
    
    grow_capacity(capacity);
    
    nodes.resize(2);
    nodes[0] = node();  // zero atom
    nodes[1] = node();  // one atom 

    /*
    0 and 1 are atoms
    0_1 and 1_0 are created right away
    0^1 has 0_1 and 1_0 as children
    0 has 0^1 as a parent 
    1 has 0^1 as a parent

        0_1 ----
        /   \     \
    0 --- 1 -- 0^1
        \   /     /
        1_0 ----

    */

    // TODO: zero out everything else


    // 0_1 sequence
    // nodes[2] = node(2, node::sequence, node_ptr(nodes, 0), node_ptr(nodes, 1));

    // // 1_0 sequence
    // nodes[3] = node(3, node::sequence, node_ptr(nodes, 1), node_ptr(nodes, 0));

    // // 0^1 relationship (refers to 0_1 and 1_0 and the parent of 0 and 1)
    // nodes[4] = node(4, node::relationship, node_ptr(nodes, 2), node_ptr(nodes, 3));
    // // nodes[0].parent = node_ptr(nodes, 4);
    // // nodes[1].parent = node_ptr(nodes, 4);

    // // 0^1 relationship
    // node_index.push_back({node_ptr(nodes, 0), node_ptr(nodes, 1)});
}


void bit_seqt::node::save(bit_seqt * s, ostream & os) const {
    os.write((const char*)&count, sizeof(count));         //  .count
    os.write((const char*)&charge, sizeof(charge));       //  .charge

    auto dex = first - &s->nodes[0];
    if(dex < 0) dex = -1;
    os.write((const char*)&dex, sizeof(dex));

    dex = second - &s->nodes[0];
    if(dex < 0) dex = -1;
    os.write((const char*)&dex, sizeof(dex));

    // dex = parent - &s->nodes[0];  
    // if(dex < 0) dex = -1;
    // os.write((const char *)&dex, sizeof(dex));
}

void bit_seqt::node::load(bit_seqt * s, istream & is) {
    is.read((char*)&count, sizeof(count));         //  .count
    is.read((char*)&charge, sizeof(charge));       //  .charge

    auto dex = first.offset();
    is.read((char*)&dex, sizeof(dex));
    if(dex < 0) first.reset();
    else first = node_ptr(s->nodes, dex);

    is.read((char*)&dex, sizeof(dex));
    if(dex < 0) second.reset();
    else second = node_ptr(s->nodes, dex);

    // is.read((char *)&dex, sizeof(dex));
    // if(dex < 0) parent.reset();
    // else parent = node_ptr(s->nodes, dex);
}


void print_seqt(bit_seqt const & s, ostream & os) {
    os << setw(12) << 0 << "| atm " << setw(12) << s.nodes[0].count << "\n";
    os << setw(12) << 1 << "| atm " << setw(12) << s.nodes[1].count << "\n";

    for(int i = 2; i < s.nodes.size(); i++) {
        bit_seqt::node const & n = s.nodes[i];
        os  << setw(12) << i;

        switch(n.op) {
        case bit_seqt::node::sequence:
            os  << "| seq " << setw(12) << n.count << " : " 
                << setw(12) << sequence_stddevs(n.first->count, n.second->count, n.count) << "% ("
                << setw(12) << n.first.offset() << " ," << setw(12) << n.second.offset() << " )\n";
            break;
        case bit_seqt::node::relationship:
            os  << "| rel " << setw(12) << n.count << " : " 
                << setw(12) << sequence_stddevs(n.first->first->count, n.first->second->count, n.first->count + n.second->count) << "%  "
                << setw(12) << n.first->first.offset() << " ~" << setw(12) << n.first->second.offset() << "  \n";
            break;
        default:
            break;
        }
    }

}


int main(int ac, char ** av) {
    bit_seqt s;

    istream * is = &cin;
    stringstream line;
    fstream f;

    if(ac > 2 && string(av[ac-2]) == "-") {
        line = stringstream(string(av[ac-1]));
        is = &line;
    } else if(ac == 2) {
        f.open(av[1]);
        if(!f) {
            cerr << "error opening '" << av[1] << "'\n";
            return -1;
        }
        is = &f;
    }

    unsigned char c, mask;

    for(;;) {
        c = is->get();
        
        if(is->eof())
            break;

        mask = 0x80;
        while(mask != 0) {
            try {
                s.read(c & mask);
            } catch(logic_error e) {
                cerr << "# Exception: " << e.what() << endl;
                goto end_read_loop;
            }
            mask >>= 1;
        }
    }
end_read_loop:
    size_t bytes_read = is->tellg();

    cerr << "# read bytes: " << bytes_read << endl;

    for(auto const & pair : s.node_index) {
        cout << "index: " << pair.first.offset() << " - " << pair.second.offset() << endl;
    }

    // fstream file("seqt.bin", ios::binary | ios::trunc | ios::out);
    // s.save(file);
    // file.close();
    print_seqt(s, cout);
    cout << endl;

    // dump out the seqt
}
