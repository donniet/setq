#include <cmath>
#include <memory>
#include <functional>
#include <map>
#include <vector>
#include <queue>
#include <tuple>
#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>

#include "char32_input_iterator.hpp"

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


double sequence_probability(double count_first, double count_second, double count_sequence) {
    // use lbinomial because this number can be very large
    double unique_sequence_count_log = lbinomial(count_first + count_second + 1, count_first);

    double sequences_with_t_transitions = 0.;
    // how many sequences have trans or higher transitions from a to b
    for(; count_sequence <= min(count_first, count_second); count_sequence++) {
        // since we are adding these we use the normal binomial function
        sequences_with_t_transitions += 
            binomial(count_sequence + 1, count_first - count_sequence)
          * binomial(count_sequence + 1, count_second - count_sequence)
        ;
    }

    return exp(log(sequences_with_t_transitions) - unique_sequence_count_log);
}

struct node;
typedef shared_ptr<node> node_ptr;

bool operator<(node_ptr const & a, node_ptr const & b) {
    return a.get() < b.get();
}

struct node {
    char32_t symbol;
    size_t count;
    double likelihood;
    node_ptr first;
    node_ptr second;

    node() 
        : symbol(0), count(0), likelihood(1), first(nullptr), second(nullptr) 
    { }
    node(char32_t s) 
        : symbol(s), count(1), likelihood(1), first(nullptr), second(nullptr) 
    { }
    node(node_ptr a, node_ptr b) 
        : symbol(0), count(1), likelihood(1), first(a), second(b) 
    { }

    bool is_atom() const { return first == nullptr && second == nullptr; }

    void tally() {
        count++;
        if(is_atom()) return;

        likelihood = sequence_probability(first->count, second->count, count);
    }
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
        return a->likelihood < b->likelihood;
    }
};

struct seqt {
    map<char32_t, node_ptr> atoms_;
    map<tuple<node_ptr, node_ptr>, node_ptr> sequences_;
    priority_queue<node_ptr, vector<node_ptr>, node_likelihood> nodes_;
    double track_likelihood_threshold;

    bool should_track(node_ptr);
    bool tally_symbol(char32_t c);

    struct reader;
};

struct seqt::reader {
    friend struct seqt;

    seqt::reader & read(char32_t c);
};



int main(int ac, char ** av) {
    char32_t c;

    ifstream f;
    istream * in = &cin;

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

    char32_input_iterator it(*in, char32_input_iterator::ascii);
    char32_input_iterator end;

    for(; it != end; ++it) { 
        c = *it;

        cout << (char)c;
    }
    cout << endl;


    return 0;
}