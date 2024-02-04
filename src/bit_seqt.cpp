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

    static float pvalue(size_t first, size_t second, size_t count_sequence);

    node_ptr ptr(size_t offset);

    void save(ostream & os);
    void read(bool bit);
    void for_each_node(std::function<void(node &)>);
    bool is_novel_suggestion(node & n);
    void pack_suggestions(node & n, node_ptr output_begin);
    void pack_suggestion_index(node & n, node_pair * output_begin);

    seqt();
};

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
    size_t index;
    size_t count;

    int charge;
    int charge2;
    float pvalue;

    node_ptr first;
    node_ptr second;
    node_ptr parent;

    node_ptr suggested_next;

    node() { } // no initialization for speed
    node(size_t dex) : 
        op(atom), index(dex), count(0), charge(0), charge2(0), 
        pvalue(0.5), first(), second(), parent(), 
        suggested_next()
    { }
    node(size_t dex, node_operation o, node_ptr f, node_ptr s) :
        op(o), index(dex), count(0), charge(0), charge2(0), 
        pvalue(0.5), first(f), second(s), parent(), 
        suggested_next()
    { }
    node(node const & r) : 
        op(r.op), index(r.index), count(r.count), charge(r.charge), charge2(r.charge2),
        pvalue(r.pvalue), first(r.first), second(r.second), parent(r.parent),
        suggested_next(r.suggested_next)
    { }
    node(node && r) : 
        op(r.op), index(r.index), count(r.count), charge(r.charge), charge2(r.charge2),
        pvalue(r.pvalue), first(r.first), second(r.second), parent(r.parent),
        suggested_next(r.suggested_next)
    { 
        r.op = atom;
        r.index = 0;
        r.count = 0;
        r.charge = 0;
        r.charge2 = 0;
        r.first.reset();
        r.second.reset();
        r.parent.reset();
    }

    node & operator=(node const & r) {
        op = r.op;
        index = r.index;
        count = r.count;
        charge = r.charge;
        charge2 = r.charge;
        first = r.first;
        second = r.second;
        parent = r.parent;
        return *this;
    }
    node & operator=(node && r) {
        op = r.op; r.op = atom;
        index = r.index;    r.index = 0;
        count = r.count;    r.count = 0;
        charge = r.charge;  r.charge = 0;
        charge2 = r.charge; r.charge2 = 0;
        first = r.first;    r.first.reset();
        second = r.second;  r.second.reset();
        parent = r.parent;  r.parent.reset();
        return *this;
    }

    // node() :
    //     index(0), count(0), charge(0), charge2(0), pvalue(0.),
    //     first(nullptr), second(nullptr), parent(nullptr)
    // { }

    void save(seqt * s, ostream & os) const;
    void load(seqt * s, istream & is);
    void compare_with(node_ptr cur);
    void multiply_charge(int c) { charge *= c; }
    void rotate_charge() { charge *= -2; } // 
    void reset_charge() { charge = 0; }
    void buffer_charge() { charge2 = charge; }
    void activate();
    void activate_atom();
    void activate_sequence();
    void activate_relationship();
    void unbuffer_charge() { charge = charge2; };
    void recalculate_pvalue();
    void clear_suggestion() { suggested_next.reset(); }

    bool has_suggestion() { return !suggested_next.is_null(); }
    bool has_parent() { return !parent.is_null(); }
    void reset_suggestion() { suggested_next.reset(); }
    bool is_active() { return charge > 1; }
};

bool seqt::node_ptr::operator<(node_ptr const & r) const {
    return &nodes_->operator[](0) + offset_ < &r.nodes_->operator[](0) + offset_;
}
bool seqt::node_ptr::operator<=(node_ptr const & r) const {
    return &nodes_->operator[](0) + offset_ <= &r.nodes_->operator[](0) + offset_;
}
seqt::node & seqt::node_ptr::operator*() const {
    return nodes_->operator[](offset_);
}
seqt::node * seqt::node_ptr::operator->() const {
    return &nodes_->operator[](offset_);
}
seqt::node_ptr::operator node*() const {
    return &nodes_->operator[](0) + offset_;
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


// TODO: this function needs to propose the creation of new nodes
void seqt::node::compare_with(node_ptr cur) {
    if(charge > 0) {
        suggested_next = cur;
    }
}

void seqt::node::recalculate_pvalue() { 
    switch(op) {
    case atom:
        pvalue = 1.;
        break;
    case sequence:
        // C(a_b)
        pvalue = seqt::pvalue(first->count, second->count, count);
        break;
    case relationship:
        // add both C(a_b) + C(b_a) to get the full count of adjacent symbols
        pvalue = seqt::pvalue(first->first->count, first->second->count, first->count + second->count);
        break;
    }
};

void seqt::node::activate_atom() {
    charge2 = 0;

    if(is_active()) {
        count++;
        charge2 = 2;
    }
}

// we use the charge2 variable to avoid race conditions
void seqt::node::activate_sequence() {
    if(charge == 1 && second->is_active()) {
        count++;
        charge2 = 2; // fully activated
    } else if(charge == 0 && first->is_active()) {
        // only add half the required charge for activation
        charge2 = 1;
    } else {
        charge2 = 0; // we didn't see the next element in the sequence
    }
}

// relationships are grandparents of their elements from the elements perspective, but parents of both
// of the possible sequences a_b and b_a
void seqt::node::activate_relationship() {
    charge2 = 0;

    if(first->first->is_active()) {
        charge2++;
        count += 2;
    }

    if(first->second->is_active()) {
        charge2++;
        count += 2;
    }
}

void seqt::node::activate() {
    // take the charge from children
    // TODO: should this be all the charge
    switch(op) {
    case atom:
        activate_atom();
        break;
    case sequence:
        activate_sequence();
        break;
    case relationship:
        activate_relationship();
        break;
    }
}



struct seqt::node_pair {
    node_ptr first;
    node_ptr second;

    bool operator<(node_pair const & r) const {
        if(first < r.first) return true;
        if(r.first < first) return false;
        return second < r.second;
    }
};


seqt::node_ptr seqt::ptr(size_t offset) {
    return node_ptr(nodes, offset);
}

void seqt::read(bool bit) {
    using namespace std::placeholders;

    size_t current_size, index_size;

    // add one to the node representing the bit read
    node_ptr cur = node_ptr(nodes, bit ? 1 : 0);

    // suggest a new node sequence node
    for_each_node(bind(&node::compare_with, _1, cur));

    // consolidate all suggestions
    //   zero out our scrap space
    //   TODO: may not be needed...
    scrap.resize(nodes.size());
    ::transform(nodes.begin(), nodes.end(), scrap.begin(), bind(&seqt::is_novel_suggestion, this, _1));
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
            bind(&seqt::pack_suggestions, this, _1, node_ptr(nodes, current_size)));

        scrap.resize(nodes.size());
        
        // update our index
        index_size = node_index.size();
        node_index.resize(index_size + suggestions);

        // add the new nodes to the index
        ::for_each(
            nodes.begin(), nodes.begin() + current_size,
            bind(&seqt::pack_suggestion_index, this, _1, &node_index[0] + index_size));

        // sort the index
        ::sort(node_index.begin(), node_index.end());
    }

    // activate charge2 using the charge variables
    for_each_node(&node::activate);

    // unbuffer the new charge2 values back into charge
    for_each_node(&node::unbuffer_charge);

    // charge the selected node
    cur->charge += 2;

    // recalculate the pvalue on all the nodes
    for_each_node(&node::recalculate_pvalue);

    // reset all the suggested_nodes
    for_each_node(&node::clear_suggestion);
    
}


void seqt::save(ostream & os) {
    for(node const & n : nodes) {
        n.save(this, os);
    }
}

bool seqt::is_novel_suggestion(node & n) {
    if(!n.has_suggestion()) 
        return false;

    node_pair suggest{node_ptr(nodes, &n-&nodes[0]), n.suggested_next};

    // ensure first is less than second
    if(&n > n.suggested_next)
        swap(suggest.first, suggest.second);

    bool exists = binary_search(node_index.begin(), node_index.end(), suggest);

    return !exists;
}

void seqt::pack_suggestions(node & n, node_ptr output_begin) {
    auto index = &n - &nodes[0]; // could use n.index?
    auto offset = output_begin.offset();
    // node at index will contain a suggestion we want to remember.
    // it will be the scrap[index]-1 index;
    if(index == 0 && scrap[index] == 0)
        return; // nothing to do

    if(index > 0 && scrap[index] == scrap[index-1])
        return; // this index does not have a novel suggestion
    
    // save the suggested node at the scrap[index]-1 index
    auto sug_index = 3*(scrap[index]-1);

    node_ptr a = ptr(offset + sug_index);
    *a = node(offset + sug_index, node::sequence, ptr(index), n.suggested_next);
    a->count = 1;
    a->charge = a->charge2 = 2;
    a->recalculate_pvalue();

    node_ptr b = ptr(offset + sug_index + 1);
    *b = node(offset + sug_index + 1, node::sequence, n.suggested_next, ptr(index));
    b->count = 0;
    b->recalculate_pvalue();

    node_ptr r = ptr(offset + sug_index + 2);
    *r = node(offset + sug_index + 2, node::relationship, ptr(offset + sug_index), ptr(offset + sug_index + 1));
    r->count = n.count + n.suggested_next->count;
    r->charge = r->charge2 = 1;

    n.parent = r;
    n.suggested_next->parent = r;

    r->recalculate_pvalue();
}

void seqt::pack_suggestion_index(node & n, node_pair * output_begin) {
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
    node_pair & a = output_begin[sug_index];
    a.first = (&n < n.suggested_next) ? pn : n.suggested_next;
    a.second = (&n < n.suggested_next) ? n.suggested_next : pn;
}

seqt::seqt() : 
    next(5), nodes(next)
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
    nodes[2] = node(2, node::sequence, node_ptr(nodes, 0), node_ptr(nodes, 1));

    // 1_0 sequence
    nodes[3] = node(3, node::sequence, node_ptr(nodes, 1), node_ptr(nodes, 0));

    // 0^1 relationship (refers to 0_1 and 1_0 and the parent of 0 and 1)
    nodes[4] = node(4, node::relationship, node_ptr(nodes, 2), node_ptr(nodes, 3));
    nodes[0].parent = node_ptr(nodes, 4);
    nodes[1].parent = node_ptr(nodes, 4);

    // 0^1 relationship
    node_index.push_back({node_ptr(nodes, 0), node_ptr(nodes, 1)});
}


void seqt::node::save(seqt * s, ostream & os) const {
    os.write((const char*)&index, sizeof(index));         //  .index
    os.write((const char*)&count, sizeof(count));         //  .count
    os.write((const char*)&charge, sizeof(charge));       //  .charge
    os.write((const char*)&pvalue, sizeof(pvalue));       //  .pvalue

    auto dex = first - &s->nodes[0];
    if(dex < 0) dex = -1;
    os.write((const char*)&dex, sizeof(dex));

    dex = second - &s->nodes[0];
    if(dex < 0) dex = -1;
    os.write((const char*)&dex, sizeof(dex));

    dex = parent - &s->nodes[0];
    if(dex < 0) dex = -1;
    os.write((const char *)&dex, sizeof(dex));
}

void seqt::node::load(seqt * s, istream & is) {
    is.read((char*)&index, sizeof(index));         //  .index
    is.read((char*)&count, sizeof(count));         //  .count
    is.read((char*)&charge, sizeof(charge));       //  .charge
    is.read((char*)&pvalue, sizeof(pvalue));       //  .pvalue

    auto dex = first.offset();
    is.read((char*)&dex, sizeof(dex));
    if(dex < 0) first.reset();
    else first = node_ptr(s->nodes, dex);

    is.read((char*)&dex, sizeof(dex));
    if(dex < 0) second.reset();
    else second = node_ptr(s->nodes, dex);

    is.read((char *)&dex, sizeof(dex));
    if(dex < 0) parent.reset();
    else parent = node_ptr(s->nodes, dex);
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

    if(ac > 2 && string(av[ac-2]) == "-") {
        line = stringstream(string(av[ac-1]));
        is = &line;
    }

    for(;;) {
        char c = is->get();
        
        if(is->eof())
            break;

        char mask = 0x80;
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
