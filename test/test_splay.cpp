
#include "splay.hpp"

int main(int, char**) {
    using std::cout, std::endl, std::string;

    splay_tree<char> * t = new splay_tree<char>();

    string str("gnarlygreenghastgah");

    for(auto i = str.begin(); i != str.end(); i++) {
        t->insert(*i);
    }

    char check [] =  { 'e',  'n',   'q',   't',   'e',  'g',  'h' };
    
    for(int i = 0; i < sizeof(check)/sizeof(char); i++) {
        cout << "string: '" << str << "' contains '" << check[i] << "'? " << t->contains(check[i]) << endl;
    }

    cout << "visiting: ";
    t->visit_preorder([](char const & c) {
        cout << c;
    });
    cout << endl;

    t->print(cout) << endl;

    for(int i = sizeof(check)/sizeof(char) - 1; i >= 0; i--) {
        cout << "deleting: " << check[i] << endl;
        t->remove(check[i]);

        t->print(cout) << endl;
    }
    
    cout << endl;

    cout << "deleting tree" << endl;
    delete t;

    return 0;
}
