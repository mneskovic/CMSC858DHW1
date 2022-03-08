#include "sparse_array/my_sparse_array.cpp" // Includes sparse array + previous two

// Runs/Demonstrates all three implementations
void test_rank_support()
{
    // Test Basic Functionality

    bit_vector simple_b = {1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0};

    my_rank_support r(&simple_b); // OG Rank support
    rank_support_v<1> r1(&simple_b); // SDSL Rank support
    string f = "rank_support/saved_rank.txt";
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
    outFile.open("rank_support/rank_test.txt");

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
    string f = "select_support/saved_select.txt";
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
    outFile.open("select_support/select_test.txt");
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
    
    string f = "sparse_array/saved_array.txt";
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
    outFile.open("sparse_array/sparse_test.txt");

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

// Runs tests
int main()
{
    test_rank_support();
    test_select_support();
    test_sparse_array();

    return 0;
}