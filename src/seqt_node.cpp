#include <cmath>
#include <memory>
#include <functional>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <tuple>
#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>

#include "utf8_reader.hpp"

/*
syntax in comments:
a precedes b: a_b
b can substitute for a: a^b

Notes:
if we counted A a's and B b's 

let's say we counted 10 a's and 10 b's.  We'd expect N(a)P(b) or 10*0.5 = 5 ab's, with a variance of 
N(a)P(b)(1-P(b)) = 2.5 (standard deviation of ~1.58)

If we instead found 9 or 10 ab's we'd suspect that a_b is an important sequence.
If we found 0 or 1 ab's we'd suspect that a and b were related in some way

in a sequence composed of C(A) a's and C(B) b's, 

C(aa) + C(ab) + C(ba) + C(bb) = C(a) + C(b) - 1

If you had 10 a's and 10 b's but we see only 2 ab sequences what is the estimated relationship between a and b?
P(a^b).

let's say P(a) = 0.5, and P(b) = 0.5, and the relationship between ab is 0.5 -- P(a^b) = 0.5

then what is the P(ab)?  P(a) * (1-P(a^b)) * (P(b) * (1-P(a^b)) + P(a) * P(a^b))



*/

using namespace std;

const char * locale = "en_US.utf8";

double binomial(double n, double k) {
    double m = max(n,k);
    double l = min(n,k);
    return tgamma(m+1) / tgamma(l+1) / tgamma(m-l+1);
}

double lbinomial(double n, double k) {
    double m = max(n,k);
    double l = min(n,k);
    return lgamma(m+1) - lgamma(l+1) - lgamma(m-l+1);
}

// what is the liklihood of p_second_given_first percent of sequences when p_first probability and p_second probability
// assuming the symbols have no relation
double sequence_probability_from_p(double p_first, double p_second, double p_firstsecond) {
    double ab = p_first * p_second;
    double v = p_second * (1-ab); // np(1-p) with n == 1/p_first
    double e = p_second / p_first; // expected number of sequences first-second
    double x = 1. / p_firstsecond; // actual number of sequences first-second
    // double x = p_first / p_second_given_first // if the third parameter is relative/given to p_first

    // this looks a lot like baysian statistics:
    // P(A|B) = P(B|A) * P(A) / P(B)

    // use erf to calculate likelihood that (P(first [then] second) - p_firstsecond) would be equal to or greater with chance
    // P-value
    return 0.5 * erfc((x-e) * sqrt(2 / v));
}

// assuming A and B are the only symbols, what is the likelihood that with 
// a = count(A) and b = count(B) that you would have ab or more sequences 
// of AB by chance.
// if this liklihood is lower than a threshold value (perhaps 1%) then 
double sequence_pvalue(uint64_t a, uint64_t b, uint64_t ab) {
    double total = (double)(a + b);
    double pa = (double)a / total;
    double pb = (double)b / total;
    double pab = pa * pb;
    double expected_value = total * pab;
    double variation = total * pab * (1-pab);

    double x = ((double)ab - expected_value) * sqrt(2 / variation);

    return 0.5 * erfc(x);
}



struct node;
struct seqt;
typedef shared_ptr<node> node_ptr;

bool operator<(node_ptr const & a, node_ptr const & b) {
    return a.get() < b.get();
}

struct node {
    char32_t symbol;
    uint64_t count;
    uint64_t sum; // sum of the sum of our children
    /*
    pvalue: this value should be around 0.5 for symbols which have no relationship
        (ie: their counts represeent the estimated values for random symbols).  
        If this pvalue is below a threshold value (say 0.01) then this sequence of symbols
        is highly unlikely to occur at random and likely represents a symbol itself.
        Conversely if the pvalue is above a threshold value (say 0.99) then this sequence
        is much more common than it would occur at random, and it is likely that these two 
        symbols have a relationship (and thus can occationally be substituted for one another,
        skeweing the counts)
    */
    double pvalue;

    /* 
    similarity should be the estimated "substitutionality" of first to second.  I believe this can be
    estimated from the sequence pvalue above.  If the pvalue is very high (like 0.99) then this sequence first_second
    appears far more likely than chance, suggesting that the two symbols can be replaced by one another.  
    I'd like to find a formula for this "replacability" based on the count of sequences.  For example, what is the expected 
    number of ab sequences when we know of A occurances of a and B occurances of b when we also expect a to be substituted
    for b 50% of the time?  It's different from the simple binomial distribution used for sequences.
    TODO: determine if this substitution is commutative (can it ever be true that a is sustituted for b 25% of the time, 
    but b is substituted for a 75% of the time?)  Is it skew symetric/complementary P(a^b) = 1-P(b^a)
    */
    // double similarity
    node_ptr first;
    node_ptr second;
    node_ptr parent; // parents are always relations

    node() 
        : symbol(0), count(0), sum(0), pvalue(1), first(nullptr), second(nullptr), parent(nullptr)
    { }
    node(char32_t s) 
        : symbol(s), count(1), sum(1), pvalue(1), first(nullptr), second(nullptr), parent(nullptr)
    { }
    node(node_ptr a, node_ptr b) 
        : symbol(0), count(1), sum(1), pvalue(1), first(a), second(b), parent(nullptr)
    { }

    bool is_atom() const { 
        return first == nullptr && second == nullptr; 
    }
    bool is_sequence() const {
        // relations have count == first->count + second->count
        // because they include everything, 
        // TODO: I think we should set first to the more frequent one, 
        // and the structure should splay or reshape itself to keep the most
        // "similar" nodes close to each other.
        return !is_atom() && pvalue < seqt::pvalue_threshold;
    }
    // TODO: perhaps relationships should have a higher bar compared
    // to sequences.  The pvalue calculation only uses the counts of the
    // symbols being considered, and ignores all the other symbols and their counts.
    // if we had a small alphabet (like binary) this would be fine, but with a large alphabet
    // like utf-8, many pairs could be considered relationships simply because of the number 
    // of other symbols that appear 
    bool is_relation() const {
        return !is_atom() && (1. - pvalue) < seqt::pvalue_threshold;
    }

    uint64_t tally() {
        count++;

        if(!is_atom())
            pvalue = sequence_pvalue(first->count, second->count, count);
    
        return count;
    }
};

/*
We need this to iterate over relationship structures.
I imagine these structures are a sort of heap, where the leaves are atoms
or sequences, and the branches are all relations.  You give 
the iterator a root relation and it iterates over all the relations, atoms, and sequences
under it.  The trick will be structuring this so that when adding new relations
there is always a root node where all the sequences and atoms that are related to a given 
node are under a single root
*/
struct node_iterator {
    // typedef node value_type;
    // typedef node& reference;
    // typedef node_ptr pointer;


    enum direction {
        first = 0, second = 1, up = 2
    };
    seqt * seqt_;
    set<node_ptr> visited_;
    vector<tuple<node_ptr, direction>> stack_;

    // pre-order traversal of relations such that for
    // (A ~ C) ~ B we visit {~, ~, A, C, B}
    bool next() {
        node_ptr p;
        direction d = up;

        // remove all the 'up's 
        while(d == up && !stack_.empty()) {
            tie(p, d) = stack_.back();
            stack_.pop_back();
            visited_.insert(p); 
        }

        if(d == up)
            return false;

        node_ptr c = (d == first) ? p->first : p->second;

        // go back up no matter what
        stack_.push_back({c, up}); // TODO: we may not need these for pre-order traversal...
        if(!visited_.contains(p) && p->is_relation()) {
            // go second if we are first, up if we are second
            // only go down if this is a relation
            stack_.push_back({c, second});
            stack_.push_back({c, first});
        }

        return true;
    }

    node_iterator & operator++() {
        next();
        return *this;
    }
    node_iterator operator++(int) {
        node_iterator r = *this;
        operator++();
        return r;
    }
    node_ptr operator->() const {
        if(stack_.empty()) return nullptr;
        return get<0>(stack_.back());
    }
    // this will error if if the stack is empty!
    node_ptr operator*() const {
        return get<0>(stack_.back());
    }
    bool operator==(node_iterator const & r) const {
        if(stack_.empty() && r.stack_.empty()) return true;  // both are empty
        if(stack_.empty() || r.stack_.empty()) return false; // one and only one is empty
        return stack_.back() == r.stack_.back();
    }
    bool operator!=(node_iterator const & r) const {
        return !(*this == r);
    }
    bool eof() const {
        return stack_.empty();
    }

    node_iterator(seqt * s, node_ptr c) 
        : seqt_(s)
    { 
        // go back up no matter what
        stack_.push_back({c, up});
        if(c->is_relation()) {
            // only go down if this is a relation
            stack_.push_back({c, second});
            stack_.push_back({c, first});
        }
    }

    node_iterator()
        : seqt_(nullptr)
    { }
};


struct node_count
{
    bool operator()(node_ptr const & a, node_ptr const & b) {
        // we use greater to put the highest count items at the top
        return a->count < b->count;
    }
};

struct node_likelihood
{
    bool operator()(node_ptr const & a, node_ptr const & b) {
        return a->pvalue < b->pvalue;
    }
};

struct seqt {
    friend struct node_collection;

    static const double pvalue_threshold = 0.003; // about 3 standard deviations

    node_ptr root_;
    
    map<char32_t, node_ptr> atoms_;
    map<tuple<node_ptr, node_ptr>, node_ptr> sequences_;
    map<node_ptr, node_ptr> rel_parent_;

    priority_queue<node_ptr, vector<node_ptr>, node_likelihood> nodes_;
    double track_likelihood_threshold;

    node_ptr atom(char32_t);
    node_ptr seq(node_ptr, node_ptr);
    node_ptr rel(node_ptr, node_ptr); // returns the lowest parent of the two parameters, or creates one

    bool should_track(node_ptr);

    struct reader;

private:
    node_ptr new_atom(uint32_t symbol);
    node_ptr new_seq(node_ptr first, node_ptr second);
    node_ptr new_rel(node_ptr first, node_ptr second);
};

node_ptr seqt::atom(char32_t c) {
    auto i = atoms_.find(c);
    if(i != atoms_.end())
        return i->second;
    
    auto r = new_atom(c);
    atoms_[c] = r;
    return r;
}

node_ptr seqt::seq(node_ptr first, node_ptr second) {
    auto i = sequences_.find({first, second});
    if(i != sequences_.end()) 
        return i->second;
    
    auto r = new_seq(first, second);
    sequences_[{first, second}] = r;
    return r;
}

node_ptr seqt::rel(node_ptr first, node_ptr second) {
    set<node_ptr> parents;
    parents.insert(first);
    parents.insert(second);
    node_ptr f = first;
    node_ptr s = second;

    while(f->parent != nullptr || s->parent != nullptr) {
        if(f->parent != nullptr) {
            f = f->parent;
            if(parents.contains(f))
                return f;
            parents.insert(f);
        }

        if(s->parent != nullptr) {
            s = s->parent;
            if(parents.contains(s))
                return s;
            parents.insert(s);
        }
    }

    // if we are here then s and f do not have parents. one of them must be root
    root_ = new_rel(f, s);
}

struct seqt::reader {
    friend struct seqt;

    seqt * seqt_;
    node_ptr context_;

    seqt::reader & read(char32_t c);
};

seqt::reader & seqt::reader::read(char32_t c) 
{
    // find or create the atom for this character
    auto a = seqt_->atom(c);
    a->tally();

    node_ptr new_context = a;

    for(node_iterator i(seqt_, context_), end; i != end; ++i) {
        node_ptr p = *i;

        auto q = seqt_->seq(p, a);
        q->tally();

        if(q->is_sequence() || q->is_relation()) {
            new_context = seqt_->rel(new_context, q);
            // don't tally because we never saw this as a sequence

        }
    }

    // replace all the read_nodes_ with the new set of read_nodes_;
    context_ = new_context;

    return *this;
}

int main(int ac, char ** av) {
    ios_base::sync_with_stdio(false);

    char32_t c;
    wifstream f;
    wistream * in = &wcin;

    // 3 sigma threshold for new symbols
    double new_symbol_threshold = 0.003;

    if(ac > 1) {
        f.open(av[1], std::ios::in);
        in = &f;
    }

    if(!*in) {
        cerr << "could not open file: " << av[1] << endl;
        return -1;
    }

    int line_count = 0;
    int char_count = 0;

    in->imbue(std::locale("C.UTF8"));

    while(true) {
        c = in->get();
        if(in->eof())
            break;

        wcout << (wchar_t)c;
    }
    wcout << endl;


    return 0;
}