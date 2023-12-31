#include "seqt.hpp"

symbol seqt::node::repr() const {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return symbol(0);

    return data->repr_;
}

seqt::node seqt::node::remove(symbol c) const {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return *this;

    if(data->ordinality_ == ordered) {
        // "removing" from an ordered node means advancing it to the next symbol if we expect c
        auto i = owner_->ordered_next_index_.find({id_, c});
        if(i == owner_->ordered_next_index_.end())
            return node(i->second, owner_); // found it, return the next ordered element

        return *this;
    }

    // "removing" from an unordered node means what you think it means
    return owner_->splay_remove(*this, c);
}

seqt::node seqt::node::expects(symbol c) {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) 
        return owner_->nil();
    
    if(is_nil())
        return owner_->nil();

    // atoms shouldn't expect anything
    // or should expect everything...
    // we should check for atomic before calling expects
    if(is_atom())
        owner_->nil();

    if(data->ordinality_ == unordered) {
        // if we are unordered, then we expect something that we haven't already seen
        if(owner_->find_in_unordered(*this, c))
            return owner_->nil();
       
        return owner_->splay_insert(*this, c);
    }


    // if we are unordered, then we should see if c 
    // we need to look through all the node's that have c as a repr_
    // and see if we are the prev of any of them.
    // I think we should create an unordered node for each ordered node
    // that contains ordered pairs of adjacent repr_ maybe?
    // i think that might be too much for now.  I think I'll cache 
    // the seqts as they are created/loaded and just look at the cache.
    // eventually I'd like the next_ id to reference an unordered collection
    // of ordered pairs of {symbol, following_id} sequences
    auto i = owner_->ordered_next_index_.find({id_, c});
    if(i != owner_->ordered_next_index_.end()) {
        return node(i->second, owner_);
    }
    return node(0, owner_);
}

bool seqt::node::is_atom() const {
    seqt_data * data = owner_->get_data(id_);

    if(data == nullptr) 
        return false;

    // unordered nodes can have zero left and right pointers
    return data->ordinality_ == ordered 
        && data->left_ == 0 
        && data->right_ == 0;
}

bool seqt::node::is_ordered() const {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return false;

    return !is_atom() && data->ordinality_ == ordered;
}

bool seqt::node::is_unordered() const {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return false;

    return !is_atom() && data->ordinality_ == unordered;
}

bool seqt::node::is_nil() const { return id_ == 0; }


seqt::node seqt::node::append(symbol c, ordinality if_atomic) const {
    if(is_nil())
        return owner_->create_atom(c);

    if(is_atom()) 
        switch(if_atomic) {
        case ordered:
            return owner_->create_ordered(*this, c);
        case unordered:
            return owner_->create_unordered(*this, c);
        }
    
    if(is_ordered()) 
        return owner_->create_ordered(*this, c);

    return owner_->create_unordered(*this, c);
}



bool seqt::node::operator==(node const & a) const {
    return id_ == a.id_ && owner_ == a.owner_;
}

bool seqt::seqt_less::operator()(seqt::node const & a, seqt::node const & b) const {
    if(a.owner_ < b.owner_) return true; // atoms are less than non-atoms
    if(b.owner_ < a.owner_) return false; 
    if(a.id_ < b.id_) return true;
    return false;
}

// don't call this on atoms, I think...
bool seqt::thread::advance(symbol c) {
    if(!unvisited().expects(c)) {
        return false;
    }

    // move this thread forward
    visited() = visited().append(c, ordered); // ensure that if visited is nil we create an atom
    unvisited() = unvisited().remove(c);
    return true;
}

seqt::node seqt::thread::visited() {
    return visited_;
}

seqt::node seqt::thread::unvisited() {
    return unvisited_;
}

seqt::seqt() : next_id_(1) {
    // insert the null node, why not?
    data_.insert(data_.begin(), seqt_data{ordered, 0, 0, 0});
}

seqt::node seqt::nil() {
    return node(0, this);
}

/* splay_remove
    This isn't exactly a splay.  I think of it like cutting a tree in
    a zig-zag pattern, and welding it together with new nodes as you go
    until you find what you are looking for.  I'm hoping this keeps the
    rough top-to-bottom order of the tree so that removal doesn't 
    elevate elements like a typical splay would.
 */
seqt::node seqt::splay_remove(seqt::node root, symbol c) {
    seqt::seqt_data * cur = get_data(root.id_);

    if(cur == nullptr || cur->ordinality_ != unordered)
        return root;

    ident ret_id = 0;
    seqt::seqt_data * tmp = nullptr;
    ident * insert_point = &ret_id;

    bool found = false;

    while(cur != nullptr) {
        if(cur->ordinality_ != unordered) {
            // it's not here
            break;
        }
        if(cur->repr_ == c) {
            // it is here
            found = true;
            break;
        }

        if(cur->repr_ < c) {
            tie(*insert_point, tmp) = new_seqt(cur->repr_, unordered);
            tmp->left_ = cur->left_;
            tmp->right_ = 0;
            insert_point = &tmp->right_;
            cur = get_data(cur->right_);
        } else { // if(c < cur->repr_)
            tie(*insert_point, tmp) = new_seqt(cur->repr_, unordered);
            tmp->right_ = cur->right_;
            tmp->left_ = 0;
            insert_point = &tmp->left_;
            cur = get_data(cur->left_);
        }
    }

    if(found) {
        // first go down the left side, but look for the first 
        // right pointer that is null and store a reference to it
        if(cur->left_ != 0) {
            *insert_point = cur->left_; // prioritize left
        }
        if(cur->right_ != 0) {
            while(*insert_point != 0) {
                // left is also non-nil, so we have to merge them 
                // by moving right until we get to a null node
                tmp = get_data(*insert_point);
                insert_point = &tmp->right_;
            }
            *insert_point = cur->right_;
        }

        return node(ret_id, this);
    }

    // we didn't find it, just return root
    // the new nodes created should get deleted since they are 
    // unlinked
    return root;
}

// returns an unordered node which is 
seqt::node seqt::splay_insert(seqt::node root, symbol c) {
    seqt::seqt_data * cur = get_data(root.id_);

    ident ret_id;
    seqt::seqt_data * ret;
    tie(ret_id, ret) = new_seqt(c, unordered);
    if(cur == nullptr)
        return node(ret_id, this);

    // let's just create new nodes for now.  
    // maybe we can filter them later.

    ident * left_insert = &ret->left_;
    ident * right_insert = &ret->right_;

    ident temp_id = root.id_;
    seqt::seqt_data * temp;

    while(cur != nullptr) {
        if(cur->ordinality_ != unordered) {
            // insert ordered seqts on the right
            *right_insert = temp_id;
            break;
        }
        if(cur->repr_ == c) {
            *left_insert = cur->left_;
            *right_insert = cur->right_;
            break;
        }
        
        if(cur->repr_ < c) {
            // we should move cur to the left
            tie(*left_insert, temp) = new_seqt(cur->repr_, unordered);
            temp->left_ = cur->left_;
            left_insert = &temp->right_;
            temp_id = cur->right_;
            cur = get_data(cur->right_);
        } 
        else { // if(c < cur->repr_)
            tie(*right_insert, temp) = new_seqt(cur->repr_, unordered);
            temp->right_ = cur->right_;
            right_insert = &temp->left_;
            temp_id = cur->left_;
            cur = get_data(cur->left_);
        }
    }

    return node(ret_id, this);
}

bool seqt::find_in_unordered(seqt::node root, symbol c) {
    seqt::seqt_data * cur = get_data(root.id_);

    while(cur != nullptr && cur->ordinality_ == unordered) {
        if(cur->repr_ < c)
            cur = get_data(cur->right_);
        else if(c < cur->repr_)
            cur = get_data(cur->left_);
        else
            return true;
    }

    return false;
}

seqt::node seqt::create_ordered(seqt::node prev, symbol c) {
    ident id;
    seqt::seqt_data * temp;
    
    auto i = ordered_next_index_.find({prev.id_, c});
    if(i != ordered_next_index_.end())
        return node(i->second, this);

    seqt::seqt_data * p = get_data(prev.id_);
    tie(id, temp) = new_seqt(c, ordered);
    if(p == nullptr)
        return node(id, this);

    temp->left_ = prev.id_;
    // add to the index
    ordered_next_index_.insert({{prev.id_, c}, id});

    return node(id, this);
}

seqt::node seqt::create_unordered(seqt::node prev, symbol c) {
    node s = splay_insert(prev, c);
    if(s.is_nil()) {
        int id;
        seqt_data * data;
        tie(id, data) = new_seqt(c, unordered);
        return node(id, this);
    }
    return s;
}

seqt::node seqt::create_atom(symbol c) {
    ident id;
    seqt_data * _;
    tie(id, _) = new_seqt(c, ordered); // atoms are considered ordered
    return node(id, this);
}

seqt::seqt_data* seqt::get_data(ident id) {
    if(id >= data_.size()) return nullptr;

    return &data_.at(id);
}

tuple<seqt::ident, seqt::seqt_data*> seqt::new_seqt(symbol c, ordinality ord) {
    ident id;
    seqt_data * data;

    auto i = atom_index_.find(c); 
    if(i != atom_index_.end())
        return {i->second, &data_[i->second]};
    
    if(!recycled_.empty()) {
        id = recycled_.back();
        recycled_.pop_back();
        data = get_data(id);

        data->ordinality_ = ord;
        data->left_ = 0;
        data->right_ = 0;
        data->repr_ = c;
    } else {
        id = next_id_++;
        data_[id] = seqt_data{ord, c, 0, 0};
        data = &data_[id];
    }

    // atom_index_.insert({c, id});
    atom_index_[c] = id;
    return {id, data};
}

void seqt::disregard(seqt::thread t) {
    // TODO: implement me!
}

void seqt::reinforce(seqt::thread t) {
    // TODO: implement me!

    threads_.push_back(t);
}


void seqt::process_threads(std::function<void(thread&)> func) {
    //HACK: this is just to avoid all the hard work of modifying the threads as we iterate over them
    thread_container temp = threads_; // copy the threads and iterate over the copy
    threads_.clear(); // clear the threads so that they are ready to be added back
    for(auto i = temp.begin(); i != temp.end(); i++)
        func(*i);
}

void seqt::read(uint32_t in) {
    auto c = symbol(in);

    process_threads([c, this](thread & t) {
        if(t.unvisited().is_atom()) {
            // if we are atomic, we only care if our representations match
            if(t.unvisited().repr() == c) {
                reinforce(t); // reinforce the atom!
            } 
        }
        else if(t.advance(c)) {
            // if we advance this thread, reinforce it
            reinforce(t);
        } 
        else {
            // otherwise disregard
            disregard(t);                                // disregard this since it didn't expect this symbol
            reinforce(thread(t.visited().append(c, ordered)));    // but start a new thread by appending the found symbol 
                                                         // onto the visited elements
        }
    });
    reinforce(thread(create_atom(c), nil()));
}
