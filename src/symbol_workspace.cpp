#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <set>
#include <functional>
#include <memory>

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
    if(c >= 33 && c <= 126) {
        os << (char)c;
        return;
    }
    os << "\\u" << std::hex << setfill('0') << setw(4) << (uint16_t)c << std::dec;
}

void escape(ostream & os, string const & s) {
    for(auto c : s) escape(os, c);
}

void write_symbol_csv(ostream & os, symbol_map const & symbols, transition_map const & transitions) {
    for(auto p : symbols) {
        os << "\"";
        escape(os, p.first);
        os << "\", " << p.second;

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


int main(int ac, char ** av) {
    istream * in = &cin;
    ifstream f;
    symbol_map counts;

    if(ac > 1) {
        f.open(av[1], std::ios::in);
        in = &f;
    }

    if(!*in) {
        cerr << "could not open file: " << av[1] << endl;
        return -1;
    }

    // count_initial_symbols(*in, counts);
    // in->seekg(std::ios::beg);

    map<string, tuple<shared_ptr<vector<string>>, int>> tracked_symbols;
    shared_ptr<vector<string>> completed = make_shared<vector<string>>();
    transition_map transitions;

    buffered_read(*in, [&](char c) {
        // increment counts on all completed symbols
        auto previously_completed = completed;
        completed = make_shared<vector<string>>();
        char sym[2] = {c, 0};
        completed->push_back(sym);

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
                for(auto prev : *prev_comp)
                    transitions[{prev, j->first}]++;

                // stop tracking this completed symbol
                tracked_symbols.erase(j);
            } else if(sym[pos] == c) {
                pos++;
            } else {
                // we didn't see the expected next symbol
                tracked_symbols.erase(j);
            }
        }

        tracked_symbols.insert({sym, {previously_completed, 1}});

        for(auto s : *completed) {
            counts[s]++;
        }
    });

    write_symbol_csv(cout, counts, transitions);

    return 0;
}
