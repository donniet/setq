#include <iostream>
#include <functional>
#include <vector>
#include <tuple>

#include <string>

/*
writing this alrogithm down with tests so I make sure I get it right
*/

template<typename T>
struct splay_tree {
    struct node {
        T value_;
        node * left_;
        node * right_;
    };

    node * root_;

    void insert(T const & v) {
        node * left = nullptr, * right = nullptr;
        node ** left_insert = &left, ** right_insert = &right;

        node * cur = root_;

        while(cur != nullptr) {
            if(cur->value_ < v) {                       // is the current node's value less than v?
                *left_insert = cur;                     //   insert our current pointer on the left
                cur = cur->right_;                      //   move our current position to the right child
                left_insert = &(*left_insert)->right_;  //   move our left insert position to the right child
                *left_insert = nullptr;                 //   unlink the right child (tracked by cur now)
            } else if(v < cur->value_) {                // otherwise do the mirror
                *right_insert = cur;                    
                cur = cur->left_;
                right_insert = &(*right_insert)->left_;
                *right_insert = nullptr;
            } else {                                    // this means we found our value
                *left_insert = cur->left_;              // insert cur's left into the left insert
                *right_insert = cur->right_;            // do the same for the right
                cur->left_ = left;                      // link our left tree with cur's left
                cur->right_ = right;                    // and our right tree with cur's right
                root_ = cur;                            // set root to cur
                return;     
            }
        }

        root_ = new node{v, left, right};
    }

    bool contains(T const & v) const {
        // we won't splay on contains to simplify things
        node * cur = root_;
        while(cur != nullptr) {
            if(cur->value_ < v)
                cur = cur->right_;
            else if(v < cur->value_)
                cur = cur->left_;
            else
                return true;
        }

        return false;
    }

    splay_tree() : root_(nullptr) { };
    ~splay_tree() {
        if(root_ == nullptr) 
            return;

        std::vector<node*> stack;
        stack.push_back(root_);

        while(!stack.empty()) {
            node * n = stack.back();
            stack.pop_back();

            if(n->left_ != nullptr) {
                stack.push_back(n->left_);
                n->left_ = nullptr;
            }
            if(n->right_ != nullptr) {
                stack.push_back(n->right_);
                n->right_ = nullptr;
            }

            delete n;
        }
    }

    void visit_preorder(std::function<void(T const &)> f) const {
        if(root_ == nullptr) return;

        std::vector<node*> stack;
        stack.push_back(root_);

        while(!stack.empty()) {
            node * n = stack.back();
            stack.pop_back();

            if(n->left_ != nullptr) stack.push_back(n->left_);
            if(n->right_ != nullptr) stack.push_back(n->right_);

            f(n->value_);
        }
    }

    void print(std::ostream & os) {
        using std::vector, std::tuple, std::tie;
        
        typedef enum {
            value = 0, comma = 1, up = 2
        } direction;

        vector<tuple<node *,direction>> stack;
        stack.push_back({root_,value});

        node * n;
        direction dir;

        while(!stack.empty()) {
            tie(n, dir) = stack.back();
            stack.pop_back();

            if(n == nullptr) {
                os << "_";
                continue;
            }

            switch(dir) {
            case value:
                os << n->value_ << "(";
                stack.push_back({n, up});
                stack.push_back({n->right_, value});
                stack.push_back({n, comma});
                stack.push_back({n->left_, value});
                break;
            case comma:
                os << ",";
                break;
            case up:
                os << ")";
                break;
            }
        }
    }
};

int main(int, char**) {
    using std::cout, std::endl, std::string;

    splay_tree<char> t;

    string str("gnarlygreenghastgah");

    for(auto i = str.begin(); i != str.end(); i++) {
        t.insert(*i);
    }

    char check [] =  { 'e',  'n',   'q',   't',   'e',  'g',  'h' };
    
    for(int i = 0; i < sizeof(check)/sizeof(char); i++) {
        cout << "string: '" << str << "' contains '" << check[i] << "'? " << t.contains(check[i]) << endl;
    }

    cout << "visiting: ";
    t.visit_preorder([](char const & c) {
        cout << c;
    });
    cout << endl;

    t.print(cout);
    cout << std::endl;

    return 0;
}
