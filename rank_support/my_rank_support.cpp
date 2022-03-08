#include <iostream>
#include <bit>
#include <cstdint>
#include "sdsl/bit_vectors.hpp" // Using sdsl-lite for bit vector
#include "sdsl/rank_support.hpp"
#include <chrono>

using namespace std;
using namespace sdsl;

class my_rank_support
{
    public:
    bit_vector * b;
    vector<uint16_t> rb;
    vector<uint32_t> rs;
    uint64_t sblock_size;
    uint64_t block_size;

    my_rank_support(bit_vector * b)
    {
        uint64_t n = b->size();

        // Get block sizes (make block_size factor of sblock_size)
        uint64_t block_size = (uint64_t)log2(n);
        this->block_size = block_size;

        uint64_t sblock_size = block_size * block_size;
        this->sblock_size = sblock_size;

        // If superblocks can't divide, resize
        int ct = n;
        while(ct % sblock_size != 0)
        {
            ct++;
        }
        b->resize(ct);

        // Set b in data structure
        this->b = b;

        n = ct;

        // Divide into superblocks
        uint64_t num_sblocks = n / sblock_size;

        vector<uint32_t> rs(num_sblocks);

        uint64_t curr_rank = 0, curr_spot = 0;
        for(uint64_t i = 0; i < num_sblocks; i++)
        {   
            for(uint64_t j = curr_spot; j < i * sblock_size; j++)
            {
                if((*b)[j] == 1)
                {
                    curr_rank += 1;
                }
            }
            rs[i] = curr_rank;
            curr_spot = i * sblock_size;
        }

        // Divide into blocks (create Rb)
        
        uint64_t num_blocks = n/block_size;

        vector<uint16_t> rb(num_blocks);
        curr_rank = 0;
        curr_spot = 0;
        for(uint64_t i = 0; i < num_blocks; i++)
        {
            for(uint64_t j = curr_spot; j < i * block_size; j++)
            {
                if((*b)[j] == 1)
                {
                    curr_rank += 1;
                }
            }
            curr_spot = i * block_size;
            rb[i] = curr_rank - rs[curr_spot / sblock_size];
        }


        // Set rs and rb
        this->rs = rs;
        this->rb = rb;

    }

    uint64_t rank1(uint64_t i)
    {
        //Using i-1 for everything to account for array vs 1 indexing

        // Get block using integer representation of binary sequence (Backwards so leftshift)
        uint64_t start = (uint64_t)((i-1)/block_size) * block_size;
        uint64_t bin_block = b->get_int(start, block_size);
        uint64_t shift = block_size - ((i-1) - start);
        uint64_t shifted_b = bin_block << ((64 - block_size) + (shift - 1));

        // Get popcount on the fly in O(1) time
        uint64_t count = popcount(shifted_b);
        uint64_t rank = (rs)[(i-1)/sblock_size] + (rb)[(i-1)/block_size] + count;    
        return rank;
    }

    uint64_t overhead()
    {
        /* Overhead is simply the size of our tables * number of bits in each item
            since popcount is on the fly
        */ 

        return 32 * rs.size() + 16 * rb.size();
    }

    void save(string &fname)
    {
        ofstream outFile;
        outFile.open(fname);

        // Save Rs
        uint64_t rs_size = rs.size();
        outFile.write((const char*)&rs_size, sizeof(rs_size));
        outFile.write((const char*)rs.data(), rs.size() * sizeof(uint32_t));

        //Save Rb
        uint64_t rb_size = rb.size();
        outFile.write((const char*)&rb_size, sizeof(rb_size));
        outFile.write((const char*)rb.data(), rb.size() * sizeof(uint16_t));

        // Save bit vector
        b->serialize(outFile);
        outFile.close();
    }

    void load(string &fname)
    {
        ifstream inFile;
        inFile.open(fname);

        // Load Rs
        uint64_t rs_size;
        rs.clear();
        inFile.read((char*)&rs_size, sizeof(rs_size));
        rs.resize(rs_size);
        inFile.read((char*)rs.data(), rs.size() * sizeof(uint32_t));

        // Load Rb
        uint64_t rb_size;
        rb.clear();
        inFile.read((char*)&rb_size, sizeof(rb_size));
        rb.resize(rb_size);
        inFile.read((char*)rb.data(), rb.size() * sizeof(uint16_t));
        
        // Load Bit Vector
        b->load(inFile);

        // Reset block sizes
        block_size = (uint64_t)log2(b->size());
        sblock_size = block_size * block_size;

        inFile.close();
    }

};