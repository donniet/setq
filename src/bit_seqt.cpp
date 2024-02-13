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

#if __clang__ && __APPLE__
// disable parallel exeuction for now
#else
#  define EXEUCUTION_PARALLEL
#endif

#if EXECUTION_PARALLEL && !DEBUG
#  define PAR_UNSEQ std::exeuction::par_unseq,
#else
#  define PAR_UNSEQ
#endif

const float statistical_significance = 3.; // standard deviations

using namespace std;

struct seqt {
    struct node;
    struct node_ptr;
    struct node_pair;

    size_t next;
    vector<node> nodes;
    vector<node_pair> node_index;
    vector<node_pair> suggested;
    vector<size_t> scrap_nodes;
    vector<size_t> scrap_index;

    vector<float> charge_buffer;

    static float pvalue(size_t first, size_t second, size_t count_sequence);

    node_ptr ptr(size_t offset);
    node_ptr ptr(node const & n);

    void save(ostream & os);
    void read(bool bit);
    void prune(size_t max_nodes);
    bool is_novel_suggestion(node_ptr first, node_ptr second);
    void pack_suggestions(node & n, node_ptr next, vector<node>::iterator output_begin);
    void pack_suggestion_index(node & n, node_ptr next, vector<node_pair>::iterator output_begin);

    seqt();
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

float seqt::pvalue(size_t first, size_t second, size_t sequence) {
    return sequence_pvalue(first, second, sequence);
}


struct seqt::node_ptr {
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


struct seqt::node {
    enum node_operation {
        atom, sequence, relationship
    } op;
    size_t count;
    int charge;

    node_ptr first;
    node_ptr second;

    node() { } // no initialization for speed
    node(node_operation o, node_ptr f, node_ptr s) :
        op(o), count(0), charge(0), first(f), second(s)
    { }
    node(node const & r) : 
        op(r.op), count(r.count), charge(r.charge), first(r.first), second(r.second)
    { }
    node(node && r) : 
        op(r.op), count(r.count), charge(r.charge), first(r.first), second(r.second)
    { 
        r.op = atom;
        r.count = 0;
        r.charge = 0;
        r.first.reset();
        r.second.reset();
    }

    node & operator=(node const & r) {
        op = r.op;
        count = r.count;
        charge = r.charge;
        first = r.first;
        second = r.second;
        return *this;
    }
    node & operator=(node && r) {
        op = r.op; r.op = atom;
        count = r.count;    r.count = 0;
        charge = r.charge;  r.charge = 0;
        first = r.first;    r.first.reset();
        second = r.second;  r.second.reset();
        return *this;
    }

    // node() :
    //     index(0), count(0), charge(0), charge2(0), pvalue(0.),
    //     first(nullptr), second(nullptr), parent(nullptr)
    // { }

    void save(seqt * s, ostream & os) const;
    void load(seqt * s, istream & is);
    // void multiply_charge(int c) { charge *= c; }
    // void rotate_charge() { charge *= -2; } // 
    // void reset_charge() { charge = 0; }
    // void buffer_charge() { charge2 = charge; }
    float significance() const;
    float activate();
    float activate_atom();
    float activate_sequence();
    float activate_relationship();
    // void unbuffer_charge() { charge = charge2; };
    // void recalculate_pvalue();
    bool is_active() { return charge > 1; }
};

bool seqt::node_ptr::operator<(node_ptr const & r) const {
    // return &nodes_->operator[](0) + offset_ < &r.nodes_->operator[](0) + r.offset_;
    // it's fine to just use the nodes_ pointer as the reference for a less than comparison
    // and it avoids a null check
    return nodes_ + offset_ < r.nodes_ + r.offset_;
    
}
bool seqt::node_ptr::operator<=(node_ptr const & r) const {
    // return &nodes_->operator[](0) + offset_ <= &r.nodes_->operator[](0) + r.offset_;
    // it's fine to just use the nodes_ pointer as the reference for a less than comparison
    // and it avoids a null check
    return nodes_ + offset_ <= r.nodes_ + r.offset_;
}
seqt::node & seqt::node_ptr::operator*() const {
    return nodes_->operator[](offset_);
}
seqt::node * seqt::node_ptr::operator->() const {
    return &nodes_->operator[](offset_);
}
seqt::node_ptr::operator node*() const {
    // we have to do a null check here I believe
    return is_null() ? nullptr : &nodes_->operator[](0) + offset_;
}
bool seqt::node_ptr::is_null() const {
    return nodes_ == nullptr || offset_ >= nodes_->size();
}
seqt::node_ptr & seqt::node_ptr::operator++() {
    if(offset_ < nodes_->size()) 
        offset_++;

    return *this;
}
seqt::node_ptr & seqt::node_ptr::operator--() {
    if(offset_ > 0)
        offset_--;

    return *this;
}
seqt::node_ptr seqt::node_ptr::operator++(int) {
    node_ptr r = *this;
    operator++();
    return r;
}
seqt::node_ptr seqt::node_ptr::operator--(int) {
    node_ptr r = *this;
    operator--();
    return r;
}

float seqt::node::significance() const {
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

float seqt::node::activate_atom() {
    if(charge > 0)
        count++;

    return 0;
}

// we use the charge2 variable to avoid race conditions
float seqt::node::activate_sequence() {
    if(charge == 1 && second->is_active()) {
        count++;
        return 2; // fully activated
    } else if(charge == 0 && first->is_active()) {
        // only add half the required charge for activation
        return 1;
    } else {
        return 0; // we didn't see the next element in the sequence
    }
}

// relationships are grandparents of their elements from the elements perspective, but parents of both
// of the possible sequences a_b and b_a
float seqt::node::activate_relationship() {
    float ret = 0;
    if(first->first->is_active()) {
        count++;
        ret = 2;
    }

    if(first->second->is_active()) {
        count++;
        ret = 2;
    }

    return ret;
}

float seqt::node::activate() {
    // take the charge from children
    // TODO: should this be all the charge
    switch(op) {
    case atom:
        return activate_atom();
    case sequence:
        return activate_sequence();
    case relationship:
        return activate_relationship();
    }
    return 0.;
}



struct seqt::node_pair {
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


seqt::node_ptr seqt::ptr(size_t offset) {
    return node_ptr(nodes, offset);
}
seqt::node_ptr seqt::ptr(node const & n) {
    return node_ptr(nodes, &n - &nodes[0]);
}

void seqt::prune(size_t max_nodes) {
    // first calculate a usefulness of each node
    vector<float> use(nodes.size());
    scrap_nodes.resize(nodes.size());

    //   get an iota equal to the number of nodes
    for_each(PAR_UNSEQ nodes.begin(), nodes.end(), [&](node const & n) {
        size_t dex = &n - &nodes[0];

        // count the dependencies;
        scrap_nodes[dex] = dex;
    });

    //   save the usefulness of each node
    for_each(PAR_UNSEQ nodes.begin(), nodes.end(), [&](node const & n) {
        size_t dex = &n - &nodes[0];

        switch(n.op) {
        case node::atom:
            use[dex] = numeric_limits<float>::max(); // infinitely useful
            break;
        case node::sequence:
            use[dex] = sequence_stddevs(n.first->count, n.second->count, n.count);
            break;
        case node::relationship:
            use[dex] = -sequence_stddevs(n.first->first->count, n.first->second->count, 
                                         n.first->count + n.second->count);
            break;
        }
    });

    // flow the usefulness values down through child nodes
    // (if a parent node is useful, then the child's usefulness should be greater)

    // sort the usefulness values and find the cutoff

    // remove the least useful nodes

    // could we combine or restructure the nodes also?  so far this code creates essentially
    // a linked list of bits.  
}

void seqt::read(bool bit) {
    using namespace std::placeholders;

    size_t current_size, index_size;

    // activate charge2 using the charge variables
    charge_buffer.resize(nodes.size());
    transform(PAR_UNSEQ nodes.begin(), nodes.end(), charge_buffer.begin(), bind(&seqt::node::activate, _1));

    // add one to the node representing the bit read
    node_ptr cur = ptr(bit ? 1 : 0);

    // suggest a new node sequence node
    scrap_nodes.resize(nodes.size());
    scrap_index.resize(nodes.size());

    for_each(PAR_UNSEQ nodes.begin(), nodes.end(), [&](node const & n) {
        size_t dex = &n - &nodes[0];
    
        if(!is_novel_suggestion(ptr(n), cur)) {
            scrap_nodes[dex] = 0;
            scrap_index[dex] = 0;
            return;
        }

        scrap_index[dex] = 1;

        if(&n == &*cur) {
            scrap_nodes[dex] = 1;
            return;
        }

        scrap_nodes[dex] = 3;
    });

    // consolidate all suggestions
    //   zero out our scrap space
    //   TODO: may not be needed...
    inclusive_scan(PAR_UNSEQ scrap_nodes.begin(), scrap_nodes.end(), scrap_nodes.begin(), plus<size_t>());
    inclusive_scan(PAR_UNSEQ scrap_index.begin(), scrap_index.end(), scrap_index.begin(), plus<size_t>());

    //   mark all the nodes with suggestions
    // ::transform_inclusive_scan(
    //     nodes.begin(), nodes.end(),
    //     scrap.begin(), 
    //     plus<size_t>(),
    //     [&](node & n) -> size_t { return is_novel_suggestion(n) ? 1 : 0; });

    // get the total suggestions from the last element of our scrap
    size_t nodes_increase = scrap_nodes.back();
    size_t index_increase = scrap_index.back();

    // now find how many individual sequences we need to add
    // inclusive_scan(PAR_UNSEQ 
    //     scrap.begin(), scrap.end(), scrap.begin(), 0, [&](size_t const & a, size_t const & b) -> size_t {
    //     return 0;
    // });


    if(nodes_increase > 0) {
        // create space for all the suggestions
        current_size = nodes.size();

        // resize the nodes to make room for the new ones
        nodes.resize(nodes.size() + nodes_increase);
        
        // pack the suggested sequences into the collection of nodes
        for_each(PAR_UNSEQ
            nodes.begin(), nodes.begin() + current_size,
            bind(&seqt::pack_suggestions, this, _1, cur, nodes.begin() + current_size));
        
        // update our index
        index_size = node_index.size();
        node_index.resize(index_size + index_increase);

        // add the new nodes to the index
        for_each(PAR_UNSEQ
            nodes.begin(), nodes.begin() + current_size,
            bind(&seqt::pack_suggestion_index, this, _1, cur, node_index.begin() + index_size));

        // sort the index
        ::sort(node_index.begin(), node_index.end());
    }

    // unbuffer the new charge2 values back into charge
    for_each(PAR_UNSEQ charge_buffer.begin(), charge_buffer.end(), [&](float const & charge) {
        size_t index = &charge - &charge_buffer[0];
        nodes[index].charge = charge;
    });

    // charge the selected node
    cur->charge += 2;
}


void seqt::save(ostream & os) {
    for(node const & n : nodes) {
        n.save(this, os);
    }
}

bool seqt::is_novel_suggestion(node_ptr first, node_ptr second) {
    if(!first->is_active() || abs(first->significance()) < statistical_significance)
        return false;

    node_pair suggest{first, second, node_ptr()};

    // if(suggest.second < suggest.first)
    //     swap(suggest.first, suggest.second);
    
    bool exists = binary_search(node_index.begin(), node_index.end(), suggest);

    return !exists;
}

void seqt::pack_suggestions(node & n, node_ptr next, vector<node>::iterator offset) {
    auto index = &n - &nodes[0]; // could use n.index?
    // node at index will contain a suggestion we want to remember.
    // it will be the scrap[index]-1 index;
    if(index == 0 && scrap_nodes[index] == 0)
        return; // nothing to do

    if(index > 0 && scrap_nodes[index] == scrap_nodes[index-1])
        return; // this index does not have a novel suggestion
    
    // save the suggested node at the scrap[index]-1 index
    int num = 1;
    if(index == 0)
        num = scrap_nodes[index];
    else
        num = scrap_nodes[index] - scrap_nodes[index-1];

    auto sug_index = scrap_nodes[index] - num;

    auto a = offset + sug_index;
    *a = node(node::sequence, ptr(n), next);
    a->count = 1;
    a->charge = 0;

    if(num == 3) {
        auto b = a + 1;
        *b = node(node::sequence, next, ptr(n));
        b->count = 0;
        // b->recalculate_pvalue();

        auto r = a + 2;
        *r = node(node::relationship, ptr(*a), ptr(*b));
        r->count = n.count + next->count;
        r->charge = 0;
    
        // n.parent = r;
        // next->parent = r;
    }

    // r->recalculate_pvalue();
}

void seqt::pack_suggestion_index(node & n, node_ptr next, vector<node_pair>::iterator output) {
    node_ptr pn(nodes, &n - &nodes[0]);
    size_t index = pn.offset();
    // node at index will contain a suggestion we want to remember.
    // it will be the scrap[index]-1 index;
    if(index == 0 && scrap_index[index] == 0)
        return; // nothing to do

    if(index > 0 && scrap_index[index] == scrap_index[index-1])
        return; // this index does not have a novel suggestion
    
    // save the suggested node at the scrap[index]-1 index
    auto sug_index = scrap_index[index]-1;

    // add this pair to the index
    auto a = output + sug_index;
    *a = node_pair(ptr(n), next, node_ptr::null());
}

seqt::seqt()
{
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
    nodes.resize(2);
    nodes[0] = node(node::atom, node_ptr::null(), node_ptr::null());  // zero atom
    nodes[1] = node(node::atom, node_ptr::null(), node_ptr::null());  // one atom 
    next = nodes.size();

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


void seqt::node::save(seqt * s, ostream & os) const {
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

void seqt::node::load(seqt * s, istream & is) {
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


void print_seqt(seqt const & s, ostream & os) {
    os << setw(12) << 0 << "| atm " << setw(12) << s.nodes[0].count << "\n";
    os << setw(12) << 1 << "| atm " << setw(12) << s.nodes[1].count << "\n";

    for(int i = 2; i < s.nodes.size(); i++) {
        seqt::node const & n = s.nodes[i];
        os  << setw(12) << i;

        switch(n.op) {
        case seqt::node::sequence:
            os  << "| seq " << setw(12) << n.count << " : " 
                << setw(12) << sequence_stddevs(n.first->count, n.second->count, n.count) << "% ("
                << setw(12) << n.first.offset() << " ," << setw(12) << n.second.offset() << " )\n";
            break;
        case seqt::node::relationship:
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
    seqt s;

    istream * is = &cin;
    stringstream line;

    auto p0 = seqt::node_pair(s.ptr(0), s.ptr(1), seqt::node_ptr::null());
    auto p2 = seqt::node_pair(s.ptr(0), s.ptr(2), seqt::node_ptr::null());
    auto p1 = seqt::node_pair(s.ptr(1), s.ptr(0), seqt::node_ptr::null());
    auto p3 = seqt::node_pair(s.ptr(1), s.ptr(2), seqt::node_ptr::null());

    if(p0.first != p1.first)
        throw std::logic_error("node_pair not ordered properly");

    if(p0 < p1 || p1 < p0)
        throw std::logic_error("node_pair less operator not working properly");

    if(p2 < p0)
        throw std::logic_error("node_pair less operator not working");

    if(p3 < p2)
        throw std::logic_error("node_pair less error");
    
    if(p3 < p0)
        throw std::logic_error("node_pair");

    if(ac > 2 && string(av[ac-2]) == "-") {
        line = stringstream(string(av[ac-1]));
        is = &line;
    }

    for(;;) {
        unsigned char c = is->get();
        
        if(is->eof())
            break;

        unsigned char mask = 0x80;
        while(mask != 0) {
            s.read(c & mask);
            mask >>= 1;
        }
    }


    // fstream file("seqt.bin", ios::binary | ios::trunc | ios::out);
    // s.save(file);
    // file.close();
    print_seqt(s, cout);
    cout << endl;

    // dump out the seqt
}
