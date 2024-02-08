#include <iostream>
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
    vector<size_t> scrap;
    vector<float> charge_buffer;

    static float pvalue(size_t first, size_t second, size_t count_sequence);

    node_ptr ptr(size_t offset);
    node_ptr ptr(node const & n);

    void save(ostream & os);
    void read(bool bit);
    void for_each_node(std::function<void(node &)>);
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
    node(size_t dex) : 
        op(atom), count(0), charge(0), first(), second()
    { }
    node(size_t dex, node_operation o, node_ptr f, node_ptr s) :
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
    return &nodes_->operator[](0) + offset_ < &r.nodes_->operator[](0) + r.offset_;
}
bool seqt::node_ptr::operator<=(node_ptr const & r) const {
    return &nodes_->operator[](0) + offset_ <= &r.nodes_->operator[](0) + r.offset_;
}
seqt::node & seqt::node_ptr::operator*() const {
    return nodes_->operator[](offset_);
}
seqt::node * seqt::node_ptr::operator->() const {
    return &nodes_->operator[](offset_);
}
seqt::node_ptr::operator node*() const {
    return nodes_ == nullptr ? nullptr : &nodes_->operator[](0) + offset_;
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
    case relationship:
        return sequence_stddevs(first->count, second->count, count);
    }
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
}



struct seqt::node_pair {
    node_ptr first;
    node_ptr second;

    node_pair() : first(), second() {}
    node_pair(node_ptr f, node_ptr s) {
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

void seqt::read(bool bit) {
    using namespace std::placeholders;

    size_t current_size, index_size;

    // activate charge2 using the charge variables
    charge_buffer.resize(nodes.size());
    ::transform(nodes.begin(), nodes.end(), charge_buffer.begin(), bind(&seqt::node::activate, _1));

    // add one to the node representing the bit read
    node_ptr cur = ptr(bit ? 1 : 0);

    // suggest a new node sequence node
    scrap.resize(nodes.size());
    ::transform(nodes.begin(), nodes.end(), scrap.begin(), [&](node const & n) {
        if(is_novel_suggestion(ptr(n), cur))
            return 1;
        return 0;
    });

    // consolidate all suggestions
    //   zero out our scrap space
    //   TODO: may not be needed...
    ::inclusive_scan(scrap.begin(), scrap.end(), scrap.begin(), plus<size_t>());

    //   mark all the nodes with suggestions
    // ::transform_inclusive_scan(
    //     nodes.begin(), nodes.end(),
    //     scrap.begin(), 
    //     plus<size_t>(),
    //     [&](node & n) -> size_t { return is_novel_suggestion(n) ? 1 : 0; });

    // get the total suggestions from the last element of our scrap
    size_t suggestions = scrap.back();

    if(suggestions > 0) {
        // create space for all the suggestions
        current_size = nodes.size();

        // resize the nodes to make room for the new ones
        nodes.resize(nodes.size() + 3*suggestions);
        
        // pack the suggested sequences into the collection of nodes
        ::for_each(
            nodes.begin(), nodes.begin() + current_size,
            bind(&seqt::pack_suggestions, this, _1, cur, nodes.begin() + current_size));
        
        // update our index
        index_size = node_index.size();
        node_index.resize(index_size + suggestions);

        // add the new nodes to the index
        ::for_each(
            nodes.begin(), nodes.begin() + current_size,
            bind(&seqt::pack_suggestion_index, this, _1, cur, node_index.begin() + index_size));

        // sort the index
        ::sort(node_index.begin(), node_index.end());
    }

    // unbuffer the new charge2 values back into charge
    ::for_each(charge_buffer.begin(), charge_buffer.end(), [&](float const & charge) {
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

    node_pair suggest{first, second};

    // if(suggest.second < suggest.first)
    //     swap(suggest.first, suggest.second);
    
    bool exists = binary_search(node_index.begin(), node_index.end(), suggest);

    return !exists;
}

void seqt::pack_suggestions(node & n, node_ptr next, vector<node>::iterator offset) {
    auto index = &n - &nodes[0]; // could use n.index?
    // node at index will contain a suggestion we want to remember.
    // it will be the scrap[index]-1 index;
    if(index == 0 && scrap[index] == 0)
        return; // nothing to do

    if(index > 0 && scrap[index] == scrap[index-1])
        return; // this index does not have a novel suggestion
    
    // save the suggested node at the scrap[index]-1 index
    auto sug_index = 3*(scrap[index]-1);

    auto a = offset + sug_index;
    *a = node(offset - nodes.begin() + sug_index, node::sequence, ptr(n), next);
    a->count = 1;
    a->charge = 2;
    // a->recalculate_pvalue();

    auto b = a + 1;
    *b = node(offset - nodes.begin() + sug_index + 1, node::sequence, next, ptr(n));
    b->count = 0;
    // b->recalculate_pvalue();

    auto r = a + 2;
    *r = node(offset - nodes.begin() + sug_index + 2, node::relationship, ptr(*a), ptr(*b));
    r->count = n.count + next->count;
    r->charge = 1;

    // n.parent = r;
    // next->parent = r;

    // r->recalculate_pvalue();
}

void seqt::pack_suggestion_index(node & n, node_ptr next, vector<node_pair>::iterator output) {
    node_ptr pn(nodes, &n - &nodes[0]);
    size_t index = pn.offset();
    // node at index will contain a suggestion we want to remember.
    // it will be the scrap[index]-1 index;
    if(index == 0 && scrap[index] == 0)
        return; // nothing to do

    if(index > 0 && scrap[index] == scrap[index-1])
        return; // this index does not have a novel suggestion
    
    // save the suggested node at the scrap[index]-1 index
    auto sug_index = scrap[index]-1;

    // add this pair to the index
    auto a = output + sug_index;
    a->first = (&n < next) ? pn : next;
    a->second = (&n < next) ? next : pn;
}

seqt::seqt() : 
    next(2), nodes(next)
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

    nodes[0] = node(0);  // zero atom
    nodes[1] = node(1);  // one atom 

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




void seqt::for_each_node(function<void(node &)> f) 
{

#if(__clang__ && __APPLE__)
    std::for_each(nodes.begin(), nodes.end(), f);
#else
    std::for_each(execution::par_unseq,
        nodes.begin(), nodes.end(),
        f);
#endif
}

template<typename OutputIterator, typename T>
void fill(OutputIterator begin, OutputIterator end, T && value)
{

#if(__clang__ && __APPLE__)
        std::fill(begin, end, value);
#else
        std::fill(execution::par_unseq, begin, end, value);
#endif

}

template<typename OutputIterator>
void sort(OutputIterator begin, OutputIterator end)
{

#if(__clang__ && __APPLE__)
        std::sort(begin, end);
#else
        std::sort(execution::par_unseq, begin, end);
#endif

}

template<typename OutputIterator, typename UnaryOp>
void for_each(OutputIterator begin, OutputIterator end, UnaryOp && value)
{

#if(__clang__ && __APPLE__)
        std::for_each(begin, end, value);
#else
        std::for_each(execution::par_unseq, begin, end, value);
#endif

}

template<typename InputIterator, typename OutputIterator, typename BinaryOp, typename UnaryOp>
void transform_inclusive_scan(InputIterator begin, InputIterator end, OutputIterator out, 
                              BinaryOp && add, UnaryOp && func)
{

#if(__clang__ && __APPLE__)
        std::transform_inclusive_scan(begin, end, out, add, func);
#else
        std::transform_inclusive_scan(execution::par_unseq, begin, end, out, add, func);
#endif

}

template<typename InputIterator, typename OutputIterator, typename UnaryOp>
void transform(InputIterator begin, InputIterator end, OutputIterator out, UnaryOp && func)
{

#if(__clang__ && __APPLE__)
        std::transform(begin, end, out, func);
#else
        std::transform(execution::par_unseq, begin, end, out, func);
#endif

}

template<typename InputIterator, typename OutputIterator, typename BinaryOp>
void inclusive_scan(InputIterator begin, InputIterator end, OutputIterator out, BinaryOp && add)
{

#if(__clang__ && __APPLE__)
        std::inclusive_scan(begin, end, out, add);
#else
        std::inclusive_scan(execution::par_unseq, begin, end, out, add);
#endif

}







int main(int ac, char ** av) {
    seqt s;

    istream * is = &cin;
    stringstream line;

    auto p0 = seqt::node_pair(s.ptr(0), s.ptr(1));
    auto p1 = seqt::node_pair(s.ptr(1), s.ptr(0));
    auto p2 = seqt::node_pair(s.ptr(0), s.ptr(2));
    auto p3 = seqt::node_pair(s.ptr(1), s.ptr(2));

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


    fstream file("seqt.bin", ios::binary | ios::trunc | ios::out);
    s.save(file);
    file.close();

    // dump out the seqt
}
