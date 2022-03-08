#include <iostream>
#include <bit>
#include <cstdint>
#include "sdsl/bit_vectors.hpp" // Using sdsl-lite for bit vector
#include "sdsl/rank_support.hpp"
#include <chrono>

using namespace std;
using namespace sdsl;

/******************Rank support class*************/
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

/*************Select Support Class**********/
class my_select_support
{

    public:
    my_rank_support * r;

    my_select_support(my_rank_support * r)
    {
        this->r = r;
    }

    uint64_t select1(uint64_t i)
    {
        uint64_t n = r->b->size();
        uint64_t start = 1;
        uint64_t end = n;
        while(start < end)
        {
            n = (start + end) / 2;
            if(r->rank1(n) > i)
            {
                end = n;
            }
            else if(r->rank1(n) == i)
            {
                if((*(r->b))[n-1] == 1)
                {
                    // If rank at position = i and is a 1, it's the ith 1
                    return n;
                }
                else
                {
                    // If it's a 0 then ith 1 comes before
                    end = n;
                }
            }
            else
            {
                start = n + 1;
            }
        }

        // Start and end overlapped which also signals ith 1
        return end;
    }

    uint64_t overhead()
    {
        // no extra space is needed, since we're just taking ranks using rank support
        return r->overhead();
    }

    // Save and load are just saving and loading rank structure
    void save(string &fname)
    {
        r->save(fname);   
    }

    void load(string &fname)
    {
        r->load(fname);   
    }
};

/**********Sparse Array Class********/
class my_sparse_array
{
    public:
    bit_vector b = bit_vector(100);
    vector<string> s;
    my_rank_support r1 = my_rank_support(&b);
    uint64_t s_size;

    void add_rank_support()
    {
        r1 = my_rank_support(&b);
    }

    void create(uint64_t size)
    {
        s_size = size;
        b.resize(size);

        // Clear b
        for(int i = 0; i < size; i++)
        {
            b[i] = 0;
        }
    }

    void append(string elem, uint64_t pos)
    {
        if(pos < b.size())
        {
            s.push_back(elem);
            b[pos] = 1;
        }
    }
    
    bool get_at_rank(uint64_t r, string & elem)
    {
        //Less than because of indexing
        if(r < s.size())
        {
            elem = s[r];
            return true;
        }
        return false;
    }

    bool get_at_index(uint64_t r, string & elem)
    {

        if(b[r] == 1)
        {
            // Using r + 1 because my rank/select use 1 indexing
            elem = s[r1.rank1(r + 1) - 1];
            b.resize(s_size);
            return true;
        }
        return false;
    }

    uint64_t num_elem_at(uint64_t r)
    {
        uint64_t num_elem = r1.rank1(r + 1);
        b.resize(s_size);
        return num_elem;
    }

    uint64_t size()
    {
        return b.size();
    }

    uint64_t num_elem()
    {
        return s.size();
    }

    void save(string & fname)
    {
        ofstream outFile;
        outFile.open(fname);

        // Save s (write each individual string on a line)
        outFile << s.size() << endl;
        for(int i = 0; i < s.size(); i++)
        {
            outFile << s[i] << endl;
        }
        // Save bit vector
        b.serialize(outFile);
        outFile.close();
    }

    void load(string & fname)
    {
        ifstream inFile;
        inFile.open(fname);

        string line;
        s.clear();

        // Load s line by line
        getline(inFile, line);
        uint64_t s_size = stoi(line);

        for(int i = 0; i < s_size; i++)
        {
            getline(inFile, line);
            s.push_back(line);
        }

        // Load bit vector
        b.load(inFile);
        inFile.close();
    }
};

void test_rank_support()
{
    // Test Basic Functionality

    bit_vector simple_b = {1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0};

    my_rank_support r(&simple_b); // OG Rank support
    rank_support_v<1> r1(&simple_b); // SDSL Rank support
    string f = "saved_rank.txt";
    r.save(f);
    bit_vector other_b = {1, 0, 1};
    my_rank_support r2(&other_b); // Loaded Rank support to test save/load
    r2.load(f);

    bool correct = true;

    for(int i = 1; i <= simple_b.size(); i++)
    {
        if(r.rank1(i) != r1.rank(i) || r.rank1(i) != r2.rank1(i))
        {
            correct = false;
        }
    }

    if(correct)
    {
        cout << "Rank works correctly" << endl;
    }
    else
    {
        cout << "Rank does not work correctly" << endl;
    }

    // Test time and space use
    ofstream outFile;
    outFile.open("rank_test.txt");

    outFile << "N,Overhead,Time" << endl;

    // Test steadily increasing N values (N squared, doubling N each iteration)
    for(uint64_t i = 4; i * i < 100000000; i *= 2){

        bit_vector b(i * i, 1); // All 1's to have varying ranks/selects
        my_rank_support r(&b);

        uint64_t rank;
        // Get time for rank operations
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rank = r.rank1((i*i) / 2);
        rank = r.rank1((i*i) / 3);
        rank = r.rank1((i*i) / 4);
        rank = r.rank1((i*i) / 5);
        rank = r.rank1((i*i) / 6);
        rank = r.rank1((i*i) / 7);
        rank = r.rank1((i*i) / 8);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        outFile << i * i << "," << r.overhead() << "," << chrono::duration_cast<chrono::nanoseconds> (end - begin).count() << endl;
    }
}

void test_select_support()
{   
    bit_vector simple_b = {1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0};

    my_rank_support r(&simple_b);
    my_select_support s(&r);
    size_t ones = rank_support_v<1>(&simple_b)(simple_b.size());
    select_support_scan<1> s1(&simple_b); // SDSL Select Support (Uses 0 indexing)
    string f = "saved_select.txt";
    s.save(f);
    bit_vector other_b = {1, 0, 1};
    my_rank_support r2(&other_b);
    my_select_support s2(&r2); // Loaded select support
    s2.load(f);
    
    bool correct = true;


    for(int i = 1; i <= 9; i++)
    {
        if(s.select1(i) != s1.select(i) + 1 || s.select1(i) != s2.select1(i))
        {
            correct = false;
        }
    }

    if(correct)
    {
        cout << "Select works correctly" << endl;
    }
    else
    {
        cout << "Select does not work correctly" << endl;
    }

    // Test time/space use
    ofstream outFile;
    outFile.open("select_test.txt");
    outFile << "N,Overhead,Time,O(logn)" << endl;
    for(int i = 4; i * i < 100000000; i *= 2){

        bit_vector b(i * i, 1); // All 1's to have varying ranks/selects
        my_rank_support r(&b);
        my_select_support s(&r);

        uint64_t select;

        // Get time for select operations
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        select = s.select1(i*i / 2);
        select = s.select1(i*i / 3);
        select = s.select1(i*i / 4);
        select = s.select1(i*i / 5);
        select = s.select1(i*i / 6);
        select = s.select1(i*i / 7);
        select = s.select1(i*i / 8);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        uint64_t time = chrono::duration_cast<chrono::nanoseconds> (end - begin).count();
        uint64_t l = (uint64_t)log2(i * i);
        outFile << i * i << "," << s.overhead() << "," << time << "," << l * 100 << endl;
    }
}

void test_sparse_array()
{
    // Basic testing
    my_sparse_array s;

    s.create(10);
    s.append("foo", 1);
    s.append("bar", 5);
    s.append("baz", 9);
    s.add_rank_support();
    
    string f = "saved_array.txt";
    s.save(f);
    my_sparse_array s1; // Loaded sparse array
    s1.load(f);

    bool correct = true;
    
    for(int i = 0; i < s.b.size(); i++)
    {
        if(s.b[i] != s1.b[i])
        {
            correct = false;
        }
    }

    for(int i = 0; i < s1.s.size(); i++)
    {
        if(s.s[i].compare(s1.s[i]) != 0)
        {
            correct = false;
        }
    }
    
    bool b_return;
    uint64_t i_return;
    string e;

    b_return = s.get_at_rank(1, e);
    if(e.compare("bar") != 0)
    {
        correct = false;
    }
    b_return = s.get_at_rank(4, e);
    if(b_return)
    {
        correct = false;
    }
    b_return = s.get_at_index(9, e);
    if(e.compare("baz") != 0)
    {
        correct = false;
    }
    b_return = s.get_at_index(6, e);
    if(b_return)
    {
        correct = false;
    }
    if(s.num_elem_at(0) != 0 || s.num_elem_at(4) != 1 || s.num_elem_at(5) != 2 || s.size() != 10 || s.num_elem() != 3)
    {
        correct = false;
    }

    if(correct)
    {
        cout << "Sparse array works correctly" << endl;
    }
    else
    {
        cout << "Sparse array does not work correctly" << endl;
    }

    ofstream outFile;
    outFile.open("sparse_test.txt");

    // Test
    
    outFile << "N,Sparsity,Time,Overhead" << endl;
    // Increase N by factor of 10
    for(int i = 100; i <= 1000000; i*=10)
    {
        my_sparse_array s2;

        // Generate sparsity
        for(int j = 1; j <= 11; j += 5)
        {
            s2.create(i);
            uint64_t wid = (uint64_t)(i / (uint64_t)(((float)j/100) * i));
            uint64_t num_1 = 0;
            for(int k = 0; k < i; k++)
            {
                if(k % wid != 0)
                {
                    num_1++;
                    s2.append("abc", k);
                }
            }

            // Needed to avoid reconstructing rank structure every time we change
            s2.add_rank_support();
            string e;

            // Get operation time
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            s2.get_at_index(i/2, e);
            s2.get_at_index(i/3, e);
            s2.get_at_rank(num_1/2, e);
            s2.get_at_rank(num_1/3, e);
            s2.num_elem_at(i/2);
            s2.num_elem_at(i/3);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            uint64_t time = chrono::duration_cast<chrono::nanoseconds> (end - begin).count();

            string es = "";
            uint64_t overhead = (i - num_1) * sizeof(es) * 8;

            outFile << i << "," << j + 1 << "," << time << "," << overhead << endl;
        }
        
    }
}
int main()
{
    test_rank_support();
    test_select_support();
    test_sparse_array();

    return 0;
}