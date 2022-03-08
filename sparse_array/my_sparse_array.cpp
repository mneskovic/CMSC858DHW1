#include "../select_support/my_select_support.cpp"

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
