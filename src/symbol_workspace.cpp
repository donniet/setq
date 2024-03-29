
// #include "partition.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <set>
#include <functional>
#include <memory>
#include <cmath>
#include <sstream>


using namespace std;


struct symbol_map 
    : public map<string, int>
{
    int & operator[](string const & s) {
        bool _;
        auto i = find(s);
        if(i == end()) 
            tie(i, _) = insert({s, 0});
            
        return i->second;
    }

    inline int & operator[](char c) {
        char s[2] = {c, 0};
        return this->operator[](s);
    }
};

struct transition_map
    : public map<pair<string,string>, int>
{
    int & operator[](pair<string,string> const & t) {
        bool _;
        auto i = find(t);
        if(i == end())
            tie(i, _) = insert({t, 0});
        
        return i->second;
    }
};


void escape(ostream & os, char c) {
    if(c == '"') {
        os << "\\\"";
    } else if(c >= 32 && c <= 126) {
        os << (char)c;
    } else {
        os << "\\u" << std::hex << setfill('0') << setw(4) << (uint16_t)c << std::dec;
    }

}

void escape(ostream & os, string const & s) {
    for(auto c : s) escape(os, c);
}

string escape_inline(string const & c) {
    stringstream ss;
    escape(ss, c);
    return ss.str();
}

void write_symbol_csv(ostream & os, symbol_map const & symbols, transition_map const & transitions) {
    os << "\"Symbol\", \"Length\", \"Count\", \"Total Characters\", \"Next...\", \"Count...\"\n";
    for(auto p : symbols) {
        os << "\"" << escape_inline(p.first) << "\", "
           << p.first.length() << ", "
           << p.second << ", "
           << p.first.length() * p.second;

        // pair<string,string> l(p.first, "");

        auto i = transitions.lower_bound({p.first, ""});
        ++i;
        vector<pair<string,int>> follows;
        for(; i->first.first == p.first; ++i) {
            follows.push_back({i->first.second, i->second});
        }

        sort(follows.begin(), follows.end(), [](auto const & a, auto const & b) {
            return a.second > b.second;
        });

        for(auto p : follows) {
            os << ", \""; escape(os, p.first); os << "\", " << p.second;
        }

        os << "\n";
    }
}

void buffered_read(istream & is, function<void(char c)> && func) {
    char buf[BUFSIZ];

    while(is) {
        is.read(buf, BUFSIZ);
        
        int num = is.gcount();
        for(int i = 0; i < num; i++) {
            func(buf[i]);
        }
    }
}

void count_initial_symbols(istream & is, symbol_map & counts) {
    buffered_read(is, [&](char c) {
        counts[c]++;
    });
}

// expects that the symbol position is not at the end of the symbol (that completed symbols have been removed)
map<string,int> read_char(char c, map<string,int> symbol_positions) {
    string s;
    int pos;

    map<string,int> ret;

    // move forward any tracked symbols
    for(auto p : symbol_positions) {
        tie(s, pos) = p;
        if(s[pos] != c)
            continue;
        
        ret[s] = pos++;
    }

    return ret;
}

map<tuple<int,int,int>, double> _probs;

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

// returns the probability that there were trans or more transitions
// in a random sequence composed of a of the first symbol and b of the second
// should be commutative in a,b
double sequence_probability(int a, int b, int trans) {
    // imagine a sequence of all a's.  we have a+1 locations to add the first b

    // baaaa, abaaa, aabaa, aaaba, aaaab                         A+1
    // bbaaaa, babaaa, baabaa, baaaba, baaaab, abbaaa, ababaa, abaaba, abaaab, aabbaa, aababa, aabaab, aaabba, aaabab, aaaabb                  
    // bbbaaaa, bbabaaa, bbaabaa, bbaaaba, babbaaa, bababaa, babaaba, babaaab,
    //          baabbaa, baababa, baabaab, baaabba, baaabab, baaaaba, 
    //          abbbaaa, abbabaa, abbaaba, abbaaab, ababbaa, abababa,
    //          ababaab, abaabba, abaabab, abaaabb, aabbbaa, aabbaba,
    //          aabbaab, aababba, aababab, aabaabb, aaabbba, aaabbab,
    //          aaababb, aaaabbb, 

    // unique_sequence_count = tgamma(ma+1) / tgamma(mi) / tgamma(ma+1 - mi);
    double unique_sequence_count_log = lbinomial(a+b+1, a);

    double sequences_with_t_transitions = 0.;
    // how many sequences have trans or higher transitions from a to b
    for(; trans <= min(a,b); trans++) {
        // start with a sequence (ab){trans} times
        // now the only places to put an a or b that won't increase the 
        // number of transitions is between an exsiting transition: aAb or aBb
        // so we have a + b - 2*trans symbols to add:
        /* abab a
           aabab abaab ababa
           abab b
           babab abbab ababb
           abab ab
           baabab babaab bababa aabbab abbaab abbaba aababb abaabb ababba 

        */

        // ____
        /* bbaaababab, baaabbabab, baaababbab, baaabababb, aaabbbabab, 
           aaabbabbab, aaabbababb, aaababbbab, aaababbabb, aaabababbb,
           bbaabaabab, bbaababaab, bbaabababa, bbabaaabab, bbabaabaab,
           bbabaababa, bbababaaab, bbababaaba, bbabababaa

        */
        // sequences_with_t_transitions += pow(trans, )
        // how many ways to put a+b things into t bins
        
        // sequences_with_t_transitions += 
        //     partition_identical_n_distinct_k(a, trans)
        //   * partition_identical_n_distinct_k(b, trans)
        // ;

        sequences_with_t_transitions += 
            binomial(trans+1, a - trans)
          * binomial(trans+1, b - trans)
        ;
    }

    return exp(log(sequences_with_t_transitions) - unique_sequence_count_log);
}



    


// getting symbols aa and bb instead of ab and ba.  Something wrong with previous completed and the transition logic
int main(int ac, char ** av) {
    istream * in = &cin;
    ifstream f;
    symbol_map counts;

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

    // count_initial_symbols(*in, counts);
    // in->seekg(std::ios::beg);

    map<string, tuple<shared_ptr<vector<string>>, int>> tracked_symbols;
    shared_ptr<vector<string>> completed = make_shared<vector<string>>();
    transition_map transitions;
    vector<transition_map::iterator>  transitions_modified;

    auto process_tracked_symbols = [&](char c = '\0') {
        // first look for transitions from the completed tracked symbols to this character
        for(auto i = tracked_symbols.begin(); i != tracked_symbols.end();) {
            auto j = i++; // save the position and increment iterator

            string const & sym = j->first;
            shared_ptr<vector<string>> prev_comp = get<0>(j->second);
            int & pos = get<1>(j->second);

            // is this symbol already completed?
            if(sym.length() <= pos) {
                // mark this symbol as completed
                completed->push_back(j->first);

                // for each previous string completed before this symbol started, count the transitions
                for(auto prev : *prev_comp) {
                    transitions[{prev, j->first}]++;
                    transitions_modified.push_back(transitions.find({prev, j->first}));
                }

                // stop tracking this completed symbol
                tracked_symbols.erase(j);
            } else if(sym[pos] == c) {
                pos++;
            } else {
                // we didn't see the expected next symbol
                tracked_symbols.erase(j);
            }
        }
    };

    auto process_transitions = [&]() {
        for(auto i : transitions_modified) {
            // cout << "transition('" << t.first.first << t.first.second << "') == " << t.second << endl;

            pair<string,string> w;
            int transition_count, a_count, b_count;

            tie(w, transition_count) = *i;    
            a_count = counts[w.first];
            b_count = counts[w.second];
            string n = w.first;
            n.append(w.second);

            // we are already tracking this transition as a symbol
            if(counts.contains(n))
                continue;

            // what is the number of ways that you can construct a sequence of a's and b's
            double transition_prob = sequence_probability(a_count, b_count, transition_count);

            if(transition_prob == 0) {
                cerr << "invalid transition probability:\n"
                     << "\t'" << escape_inline(w.first) << "': " << a_count << "\n"
                     << "\t'" << escape_inline(w.second) << "': " << b_count << "\n"
                     << "\ttransitions: " << transition_count << endl;
                throw logic_error("transition probability is zero");
            }

            // after a symbol is added we should stop tracking the transition altogether maybe?
            if(transition_prob < new_symbol_threshold && !counts.contains(n)) {
                //cout << "new symbol: " << escape_inline(n) << " p-value: " << setw(5) << transition_prob * 100. << "%" << endl;
                
                /*
                    sequence_node nn(n);
                    nn.first = w.first;
                    nn.second = w.second;
                */
                counts[n] = transition_count;
            }
            /*
                or should we check to see if this transition crosses back into the "very likely" like < 1sigma
            */
        }
        transitions_modified.clear();
    };
    
    stringstream contents;

    while(*in) {
        char c = in->get();
        contents.put(c);
        
        if(c == '\n') {
            line_count++;
            char_count = 0;
        } else {
            char_count++;
        }

        // increment counts on all completed symbols
        completed = make_shared<vector<string>>();
        char s[2] = {c, 0};

        // completed->push_back(s);

        process_tracked_symbols(c);

        for(auto s : *completed) {
            counts[s]++;
        }

        auto i = counts.lower_bound(s);
        if(i != counts.end()) i++;
        for(; i != counts.end() && i->first[0] == c; i++) {          
            // cerr << "tracking: " <<  i->first << endl;
            tracked_symbols.insert({i->first, {completed, 1}});
        }

        // don't add sym because we already marked it as completed
        // only add tracked_symbols if they aren't completed.
        tracked_symbols.insert({s, {completed, 1}});

        process_transitions();
    };

    // one more time to close out anything
    completed = make_shared<vector<string>>();
    process_tracked_symbols();
    for(auto s : *completed) {
        counts[s]++;
    }
    process_transitions();

    /* now lets figure out which are potential new symbols to create.
     these are transitions very unlikely to happen by chance
    */

    // for(auto p : counts) {
    //     cout << "count('" << p.first << "') == " << p.second << endl; 
    // }

    // run through the file again, but zero out the counts
    counts.clear();
    contents.seekg(ios::beg);
    while(contents) {
        char c = contents.get();

        // increment counts on all completed symbols
        completed = make_shared<vector<string>>();
        char s[2] = {c, 0};

        // completed->push_back(s);

        process_tracked_symbols(c);

        for(auto s : *completed) {
            counts[s]++;
        }

        auto i = counts.lower_bound(s);
        if(i != counts.end()) i++;
        for(; i != counts.end() && i->first[0] == c; i++) {          
            // cerr << "tracking: " <<  i->first << endl;
            tracked_symbols.insert({i->first, {completed, 1}});
        }

        // don't add sym because we already marked it as completed
        // only add tracked_symbols if they aren't completed.
        tracked_symbols.insert({s, {completed, 1}});

        process_transitions();
        // transitions_modified.clear();  // we won't process transitions because we don't want new symbols this time

    }

    write_symbol_csv(cout, counts, transitions);

    return 0;
}
