#include "seqt.hpp"


mem::seqt mem::seqt::remove(symbol c) const {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return *this;

    if(data->ordinality_ == ordered) {
        // "removing" from an ordered seqt means advancing it to the next symbol if we expect c
        auto i = owner_->ordered_next_index_.find({id_, c});
        if(i == owner_->ordered_next_index_.end())
            return seqt(i->second, owner_); // found it, return the next ordered element

        return *this;
    }

    // "removing" from an unordered seqt means what you think it means
    return owner_->splay_remove(*this, c);
}

mem::seqt mem::seqt::expects(symbol c) {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return false;

    if(data->ordinality_ == unordered) {
        // if we are unordered, then we expect something that we haven't already seen
        return !owner_->find_in_unordered(*this, c);
    }

    // if we are unordered, then we should see if c 
    // we need to look through all the seqt's that have c as a repr_
    // and see if we are the prev of any of them.
    // I think we should create an unordered seqt for each ordered seqt
    // that contains ordered pairs of adjacent repr_ maybe?
    // i think that might be too much for now.  I think I'll cache 
    // the seqts as they are created/loaded and just look at the cache.
    // eventually I'd like the next_ id to reference an unordered collection
    // of ordered pairs of {symbol, following_id} sequences
    auto i = owner_->ordered_next_index_.find({id_, c});
    if(i != owner_->ordered_next_index_.end()) {
        return seqt(i->second, owner_);
    }
    return seqt(0, owner_);
}

bool mem::seqt::is_ordered() const {
    seqt_data * data = owner_->get_data(id_);
    if(data == nullptr) return false;

    return data->ordinality_ == ordered;
}

mem::seqt mem::seqt::append(symbol c, bool default_ordered) const {
    if(is_nil()) {
        if(default_ordered) return owner_->create_ordered(c);
        return owner_->create_unordered(c);
    }
    if(is_ordered()) return owner_->create_ordered(*this, c);
    return owner_->create_unordered(*this, c);
}

bool mem::seqt::is_nil() const { return id_ == 0; }


bool mem::seqt::operator==(seqt const & a) const {
    return id_ == a.id_ && owner_ == a.owner_;
}

bool mem::seqt_less::operator()(mem::seqt const & a, mem::seqt const & b) const {
    if(a.owner_ < b.owner_) return true; // atoms are less than non-atoms
    if(b.owner_ < a.owner_) return false; 
    if(a.id_ < b.id_) return true;
    return false;
}

bool mem::thread::advance(symbol c) {
    if(!unvisited().expects(c)) {
        return false;
    }

    // move this thread forward
    visited() = visited().append(c);
    unvisited() = unvisited().remove(c);
    return true;
}

mem::seqt mem::thread::visited() {
    return visited_;
}

mem::seqt mem::thread::unvisited() {
    return unvisited_;
}

mem::mem() : next_id_(1) {
    // insert the null seqt, why not?
    data_.insert(data_.begin(), seqt_data{ordered, 0, 0, 0});
}

mem::seqt mem::nil() {
    return seqt(0, this);
}

/* splay_remove
    This isn't exactly a splay.  I think of it like cutting a tree in
    a zig-zag pattern, and welding it together with new nodes as you go
    until you find what you are looking for.  I'm hoping this keeps the
    rough top-to-bottom order of the tree so that removal doesn't 
    elevate elements like a typical splay would.
 */
mem::seqt mem::splay_remove(mem::seqt root, symbol c) {
    mem::seqt_data * cur = get_data(root.id_);

    if(cur == nullptr || cur->ordinality_ != unordered)
        return root;

    ident ret_id = 0;
    mem::seqt_data * tmp = nullptr;
    ident * insert_point = &ret_id;

    /*      @        @ is new
           / \  <--- * is copied
          *   x      x is the new insertion point
     */
    auto copy_left = [&](mem::seqt_data * l) { 
        tie(*insert_point, tmp) = new_seqt(l->repr_, unordered);
        tmp->left_ = l->left_;
        tmp->right_ = 0;
        insert_point = &tmp->right_;
    };

    /*      @        @ is new
           / \  <--- * is copied
          x   *      x is the new insertion point
     */
    auto copy_right = [&](mem::seqt_data * r) { 
        tie(*insert_point, tmp) = new_seqt(r->repr_, unordered);
        tmp->right_ = r->right_;
        tmp->left_ = 0;
        insert_point = &tmp->left_;    
    };
    
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

        return seqt(ret_id, this);
    }

    // we didn't find it, just return root
    // the new nodes created should get deleted since they are 
    // unlinked
    return root;
}

// returns an unordered seqt which is 
mem::seqt mem::splay_insert(mem::seqt root, symbol c) {
    mem::seqt_data * cur = get_data(root.id_);

    ident ret_id;
    mem::seqt_data * ret;
    tie(ret_id, ret) = new_seqt(c, unordered);
    if(cur == nullptr)
        return seqt(ret_id, this);

    // let's just create new nodes for now.  
    // maybe we can filter them later.

    ident * left_insert = &ret->left_;
    ident * right_insert = &ret->right_;

    ident temp_id = root.id_;
    mem::seqt_data * temp;

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

    return seqt(ret_id, this);
}

bool mem::find_in_unordered(mem::seqt root, symbol c) {
    mem::seqt_data * cur = get_data(root.id_);

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

mem::seqt mem::create_ordered(mem::seqt prev, symbol c) {
    ident id;
    mem::seqt_data * temp;
    
    auto i = ordered_next_index_.find({prev.id_, c});
    if(i != ordered_next_index_.end())
        return seqt(i->second, this);

    mem::seqt_data * p = get_data(prev.id_);
    tie(id, temp) = new_seqt(c, ordered);
    if(p == nullptr)
        return seqt(id, this);

    temp->left_ = prev.id_;
    // add to the index
    ordered_next_index_.insert({{prev.id_, c}, id});

    return seqt(id, this);
}

mem::seqt mem::create_unordered(mem::seqt prev, symbol c) {
    seqt s = splay_insert(prev, c);
    if(s.is_nil()) {
        int id;
        seqt_data * data;
        tie(id, data) = new_seqt(c, unordered);
        return seqt(id, this);
    }
    return s;
}

mem::seqt mem::create_ordered(symbol c) {
    ident id;
    seqt_data * _;
    tie(id, _) = new_seqt(c, ordered);
    return seqt(id, this);
}

mem::seqt mem::create_unordered(symbol c) {
    ident id;
    seqt_data * _;
    tie(id, _) = new_seqt(c, unordered);
    return seqt(id, this);
}

mem::seqt_data* mem::get_data(ident id) {
    if(id >= data_.size()) return nullptr;

    return &data_.at(id);
}

tuple<mem::ident, mem::seqt_data*> mem::new_seqt(symbol c, ordinality ord) {
    ident id;
    seqt_data * data;

    switch(ord) {
    case ordered:
        // use a for loop for scope
        for(auto i = ordered_atom_index_.find(c); i != ordered_atom_index_.end();)
            return {i->second, &data_[i->second]};
        break;
    case unordered:
        for(auto i = unordered_atom_index_.find(c); i != unordered_atom_index_.end();)
            return {i->second, &data_[i->second]};
        break;
    }

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
        data = &*data_.insert(data_.begin() + id, seqt_data{
            ord, c, 0, 0
        }); 
    }
    return {id, data};
}

void mem::disregard(mem::thread t) {
    // TODO: implement me!
}

void mem::reinforce(mem::thread t) {
    // TODO: implement me!

    threads_.push_back(t);
}


void mem::process_threads(std::function<void(thread&)> func) {
    //HACK: this is just to avoid all the hard work of modifying the threads as we iterate over them
    thread_container temp = threads_; // copy the threads and iterate over the copy
    threads_.clear(); // clear the threads so that they are ready to be added back
    for(auto i = temp.begin(); i != temp.end(); i++)
        func(*i);
}

void mem::read(uint32_t in) {
    auto c = symbol(in);

    process_threads([c, this](thread & t) {
        if(t.advance(c)) {
            reinforce(t);
        } 
        else {
            disregard(t);                                // disregard this since it didn't expect this symbol
            reinforce(thread(t.visited().append(c)));    // but start a new thread by appending the found symbol 
                                                         // onto the visited elements
        }
    });
    reinforce(thread(create_ordered(c), nil()));
    reinforce(thread(create_unordered(c), nil()));
}
