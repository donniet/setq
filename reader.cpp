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
    int charge;

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

    int sig() const {
        if(charge > threshold) return 1;
        if(charge < threshold) return -1;
        return 0;
    }
};

struct seqt {
    std::shared_mutex read_mutex;
    std::condition_variable_any cv;

    // these members must only be modified under a unique_lock
    std::map<symbol, identifier> symbol_index;
    identifier next_id;
    bool closed;
    size_t symbols_read;
    typedef enum {
        calculate_charges  = 0,
        pump_charges       = 1
    } stage;

    stage current_stage;

    // these two methods demand a unique_lock
    void read(symbol c);
    void close();


    std::vector<seqt_entry> data_;
    std::vector<size_t> charge_buffer_;


    void increase_buffer(uint64_t max_size);
    void initialize_seqt_entry_atom(seqt_entry &, symbol c);

    // thinking thread works on `ID % modulus == remainder` items
    // this method is designed as a worker and demands a shared_lock
    void think(identifier modulus, identifier remainder);

    int do_calculate_charges(identifier id);
};

int seqt::do_calculate_charges(identifier id) {
    auto repr = data_[id].repr;

    int ret = 0;

    if(data_[id].type == atom) {
        
    } 
    else if(data_[id].type == sequence) {
        auto seq_prev = data_[id].seq_prev;
        // is the previous sequence reseting?  Then we should reset regardless of our repr
        if(data_[seq_prev].charge < data_[seq_prev].threshold)
            charge_buffer_[id]--;
        // are we at equilibrium or firing and the previous sequence is equilibrium or firing?  Then we should fire
        else if(data_[repr].charge >= data_[repr].threshold && 
                data_[seq_prev].charge >= data_[seq_prev].threshold) 
            charge_buffer_[id]++;
    }
    else if(data_[id].type == set) {
        auto set_left = data_[id].set_left;
        auto set_right = data_[id].set_right;

        if( data_[repr].charge >= data_[repr].threshold && 
            data_[set_left].charge >= data_[set_left].threshold &&
            data_[set_right].charge >= data_[set_right].threshold)
            
            charge_buffer_[id]++;

        // maybe we shouldn't reset sets since we don't know how big they are?

    }
}

void seqt::think(identifier modulus, identifier remainder) {
    std::shared_lock<std::shared_mutex> shared_lock(read_mutex);
    size_t my_symbols_read = symbols_read;
    stage my_stage = current_stage;

    for(;;) {
        cv.wait(shared_lock, [&my_symbols_read, &my_stage, this](){
            return closed || my_symbols_read != symbols_read || my_stage != current_stage;
        }); // wait to be notified that a read occured

        if(closed) break; // we are done.

        switch(current_stage) {
        case calculate_charges:
            for(identifier id = remainder; id < next_id; id += modulus) {
                do_calculate_charges(id);
            }
            break;
        case pump_charges:
            for(identifier id = remainder; id < next_id; id += modulus) {
                data_[id].charge += charge_buffer_[id];
            }
            break;
        }

        my_symbols_read = symbols_read;
    }
}

void seqt::close() {
    std::unique_lock<std::shared_mutex> lock(read_mutex);
    closed = true;
    lock.unlock();
    cv.notify_all();
}

void seqt::read(symbol c) {
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
    data_[symbol_id].charge++;

    lock.unlock();
    cv.notify_all();
}