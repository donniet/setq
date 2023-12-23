#ifndef __SPLAY_HPP__
#define __SPLAY_HPP__

#include <iostream>
#include <functional>
#include <vector>
#include <tuple>

#include <string>

using std::vector, std::tuple, std::tie, std::cerr, std::cout, 
      std::endl, std::ostream;
/*
writing this alrogithm down with tests so I make sure I get it right
*/

template<typename T>
struct splay_tree 
{
    struct node 
    {

        T value_;
        node * left_;
        node * right_;
        size_t id_;

        node(T value, node * left, node * right, size_t id)
            : value_(value), left_(left), right_(right), id_(id)
        { 
            cerr << "create: " << id_ << endl;
        }

        ~node() 
        {
            cerr << "destroy: " << id_ << endl;
        }
    };

    node * root_;
    size_t next_id_;

    void remove(T const & v) 
    {
        node * left = nullptr, * right = nullptr;
        node ** left_insert = &left, ** right_insert = &right;

        node * cur = root_;
        root_ = nullptr; // unlink root for accounting purposes (maybe smart pointers?)

        bool found = false;

        while(cur != nullptr) 
        {
            if(cur->value_ < v) 
            {
                *left_insert = cur;
                cur = cur->right_;
                left_insert = &(*left_insert)->right_;
                *left_insert = nullptr;
            }
            else if(v < cur->value_) 
            {
                *right_insert = cur;
                cur = cur->left_;
                right_insert = &(*right_insert)->left_;
                *right_insert = nullptr;
            }
            else 
            {
                found = true;
                break;
            }
        }

        if(found) 
        {
            if(cur->left_ != nullptr) 
            {
                *left_insert = cur->left_;
                cur->left_ = nullptr;
            }
            if(cur->right_ != nullptr) 
            {
                *right_insert = cur->right_;
                cur->right_ = nullptr;
            }
            // deleting the found node
            delete cur;
        }

        root_ = merge(left, right);
    }

    static node * merge(node * left, node * right, bool prefer_left = true) 
    {
        if(left == nullptr)
            return right; // could be null

        if(right == nullptr)
            return left;
        
        // both are non-null
        if(prefer_left)
        {
            node * max = max_node(left);
            max->right_ = right;
            return left;
        }
        node * max = max_node(right);
        max->right_ = left;
        return right;

    }

    // does not check for null
    static node * max_node(node * r) 
    {
        while(r->right_ != nullptr)
            r = r->right_;

        return r;
    }

    void insert(T const & v) 
    {
        node * left = nullptr, * right = nullptr;
        node ** left_insert = &left, ** right_insert = &right;

        node * cur = root_;

        while(cur != nullptr) {
            if(cur->value_ < v) {                       // is the current node's value less than v?
                *left_insert = cur;                     //   insert our current pointer on the left
                cur = cur->right_;                      //   move our current position to the right child
                left_insert = &(*left_insert)->right_;  //   move our left insert position to the right child
                *left_insert = nullptr;                 //   unlink the right child (tracked by cur now)
            } 
            else if(v < cur->value_) {                  // otherwise do the mirror
                *right_insert = cur;                    
                cur = cur->left_;
                right_insert = &(*right_insert)->left_;
                *right_insert = nullptr;
            } 
            else {                                      // this means we found our value
                *left_insert = cur->left_;              // insert cur's left into the left insert
                *right_insert = cur->right_;            // do the same for the right
                cur->left_ = left;                      // link our left tree with cur's left
                cur->right_ = right;                    // and our right tree with cur's right
                root_ = cur;                            // set root to cur
                return;     
            }
        }

        root_ = new node(v, left, right, next_id_++);
    }

    bool contains(T const & v) const 
    {
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

    splay_tree() 
        : root_(nullptr), next_id_(0) 
    { };

    ~splay_tree()
    {
        destroy_tree(root_);
    }

    static void destroy_tree(node *& root) 
    {
        if(root == nullptr) 
            return;

        vector<node*> stack;
        stack.push_back(root);

        while(!stack.empty()) 
        {
            node * n = stack.back();
            stack.pop_back();

            if(n->left_ != nullptr) 
            {
                stack.push_back(n->left_);
                n->left_ = nullptr;
            }
            if(n->right_ != nullptr) 
            {
                stack.push_back(n->right_);
                n->right_ = nullptr;
            }

            delete n;
        }

        root = nullptr;
    }

    void visit_preorder(std::function<void(T const &)> func) const 
    {
        if(root_ == nullptr) return;

        vector<node*> stack;
        stack.push_back(root_);

        while(!stack.empty()) {
            node * n = stack.back();
            stack.pop_back();

            if(n->left_ != nullptr) 
                stack.push_back(n->left_);
            if(n->right_ != nullptr) 
                stack.push_back(n->right_);

            func(n->value_);
        }
    }

    ostream & print(ostream & os) 
    {
        if(root_ == nullptr) 
        {
            os << "\n";
            return os;
        }
        
        typedef enum 
        {
            value = 0, comma = 1, up = 2
        } direction;

        vector<tuple<node *,direction, int>> stack;
        stack.push_back({root_,value, 0});

        node * n;
        direction dir;
        int depth;

        while(!stack.empty()) 
        {
            tie(n, dir, depth) = stack.back();
            stack.pop_back();

            if(n == nullptr) 
            {
                for(int d = depth; d >= 0; d--) 
                    os << " ";
                os << "_\n";

                continue; // while
            }

            switch(dir) {
            case value:
                for(int d = depth; d >= 0; d--) 
                    os << " ";

                os << n->value_ << "\n";
                stack.push_back({n, up, depth});
                if(n->right_ != nullptr) 
                    stack.push_back({n->right_, value, depth+1});

                stack.push_back({n, comma, depth});
                if(n->left_ != nullptr)
                    stack.push_back({n->left_, value, depth+1});

                break;
            case comma:
                // os << ",";
                break;
            case up:
                // os << ")";
                break;
            }
        }
        return os;
    }
};


#endif
