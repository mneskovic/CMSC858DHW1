#include "../rank_support/my_rank_support.cpp"

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
