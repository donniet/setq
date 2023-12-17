#include <shared_mutex>
#include <condition_variable>
#include <vector>
#include <map>

typedef uint32_t symbol;
typedef uint64_t identifier;

typedef enum {
    nil = 0, atom, sequence, set
} seqt_type;

struct seqt_entry {
    seqt_type type;
    identifier id;

    size_t threshold;
    size_t weight;

    // atom
    symbol c;

    // non-atom
    identifier repr;

    // sequence
    identifier seq_prev;
    identifier seq_next;

    // tree
    identifier set_left;
    identifier set_right;
};

struct mem {
    std::shared_mutex read_mutex;
    std::condition_variable_any cv;

    // these members must only be modified under a unique_lock
    std::map<symbol, identifier> symbol_index;
    identifier next_id;
    bool closed;
    size_t symbols_read;

    // these two methods demand a unique_lock
    void read(symbol c);
    void close();


    std::vector<seqt_entry> data_;
    void increase_buffer(uint64_t max_size);
    void initialize_seqt_entry_atom(seqt_entry &, symbol c);

    // thinking thread works on `ID % modulus == remainder` items
    // this method is designed as a worker and demands a shared_lock
    void think(identifier modulus, identifier remainder);
};

void mem::think(identifier modulus, identifier remainder) {
    std::shared_lock<std::shared_mutex> lock(read_mutex);
    size_t my_symbols_read = symbols_read;

    // actually this won't work in parallel...
    std::vector<identifier> fires;
    std::vector<identifier> resets;

    for(;;) {
        cv.wait(lock, [my_symbols_read, this](){
            return closed || my_symbols_read != symbols_read;
        }); // wait to be notified that a read occured

        if(closed) break; // we are done.

        // here's where we have to decide how to spread the weights out
        // in our structure.  At first maybe we just propogate them?
        for(identifier id = remainder; id < next_id; id += modulus) {
            auto repr = data_[id].repr;

            if(data_[id].type == sequence) {
                auto seq_prev = data_[id].seq_prev;
                // is the previous sequence reseting?  Then we should reset regardless of our repr
                if(data_[seq_prev].weight < data_[seq_prev].threshold)
                    resets.push_back(id);
                // are we at equilibrium or firing and the previous sequence is equilibrium or firing?  Then we should fire
                else if(data_[repr].weight >= data_[repr].threshold && 
                        data_[seq_prev].weight >= data_[seq_prev].threshold) 
                    fires.push_back(id);
            }
            else if(data_[id].type == set) {
                auto set_left = data_[id].set_left;
                auto set_right = data_[id].set_right;

                if( data_[repr].weight >= data_[repr].threshold && 
                    data_[set_left].weight >= data_[set_left].threshold &&
                    data_[set_right].weight >= data_[set_right].threshold)
                    
                    fires.push_back(id);

                // maybe we shouldn't reset sets since we don't know how big they are?

            }

        }


        my_symbols_read = symbols_read;
        fires.clear();
    }
}

void mem::close() {
    std::unique_lock<std::shared_mutex> lock(read_mutex);
    closed = true;
    lock.unlock();
    cv.notify_all();
}

void mem::read(symbol c) {
    std::unique_lock<std::shared_mutex> lock(read_mutex);

    // get the id of this symbol
    identifier symbol_id = 0;
    auto i = symbol_index.find(c);
    if(i == symbol_index.end()) {
        if(next_id >= data_.size())
            increase_buffer(next_id);
        symbol_id = next_id++;
        symbol_index[c] = symbol_id;
        initialize_seqt_entry_atom(data_[symbol_id], c);
    } else {
        symbol_id = i->second;
    }

    // increment the weight of this symbol
    data_[symbol_id].weight++;

    lock.unlock();
    cv.notify_all();
}