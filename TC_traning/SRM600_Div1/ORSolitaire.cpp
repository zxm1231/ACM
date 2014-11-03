#line 2 "ORSolitaire.cpp"
#include <bits/stdc++.h>

using namespace std;
#define outstars cout << "***********************" << endl;
#define clr(a,b) memset(a,b,sizeof(a))
#define mk make_pair
#define pb push_back
#define sz size()
#define AA first
#define BB second
#define eps 1e-10
#define zero(x) fabs(x) < eps
#define deq(a , b) fabs(a - b) < eps

const int MAXN = 40000 + 50;
const int MAXE = 4000 + 50;
const int MAXQ = 1000000 + 50;

const int inf = 0x3f3f3f3f;
const int INF = ~0U >> 1;
const long long LLinf = 0x3FFFFFFFFFFFFFFFLL;
const long long LLINF = (1LL << 63) - 1;
const int IMIN = 0x80000000;

const int mod = 10007;
const long long MOD = 1000000000 + 7;
const double PI = acos(-1.0);

typedef long long LL;
typedef pair<int , int> pii;
typedef vector<int> vec;


int a[35];
int cnt[35];
class ORSolitaire
{
public:
    int getMinimum(vector <int> numbers, int goal)
    {
        clr(a , 0);
        clr(cnt , 0);
        bitset <32> g(goal);
        string sg = g.to_string();
        for(int i = 0 ; i < 32 ; i ++) {
            if(sg[i] == '1')a[i] = 1;
        }
        int ans = INF;
        for(int i = 0 ; i < numbers.sz ; i++) {
            bitset <32> n(numbers[i]);
            string sn = n.to_string();
            int ok = 1;
            for(int i = 0 ; i < 32 ; i++) {
                if(sn[i] == '1' && !a[i]) {
                    ok = 0;
                    break;
                }
            }
            if(!ok)continue;
            for(int i = 0 ; i < 32 ; i++) {
                if(sn[i] == '1')cnt[i]++;
            }
        }
//        for(int i = 0 ; i < 32 ; i++) {
//            cout << i << ' ' << a[i] << endl;
//        }
//        for(int i = 0 ; i < 32 ; i++) {
//            cout << i << ' ' << cnt[i] << endl;
//        }
        for(int i = 0 ; i < 32 ; i++) {
            if(a[i])ans = min(ans , cnt[i]);
        }
        return ans;
    }

// BEGIN CUT HERE
public:
    void run_test(int Case)
    {
        if ((Case == -1) || (Case == 0)) test_case_0();
        if ((Case == -1) || (Case == 1)) test_case_1();
        if ((Case == -1) || (Case == 2)) test_case_2();
        if ((Case == -1) || (Case == 3)) test_case_3();
        if ((Case == -1) || (Case == 4)) test_case_4();
    }
private:
    template <typename T> string print_array(const vector<T> &V)
    {
        ostringstream os;
        os << "{ ";
        for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\",";
        os << " }";
        return os.str();
    }
    void verify_case(int Case, const int &Expected, const int &Received)
    {
        cerr << "Test Case #" << Case << "...";
        if (Expected == Received) cerr << "PASSED" << endl;
        else {
            cerr << "FAILED" << endl;
            cerr << "\tExpected: \"" << Expected << '\"' << endl;
            cerr << "\tReceived: \"" << Received << '\"' << endl;
        }
    }
    void test_case_0()
    {
        int Arr0[] = {1, 2, 4};
        vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
        int Arg1 = 7;
        int Arg2 = 1;
        verify_case(0, Arg2, getMinimum(Arg0, Arg1));
    }
    void test_case_1()
    {
        int Arr0[] = {1, 2, 4, 7, 8};
        vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
        int Arg1 = 7;
        int Arg2 = 2;
        verify_case(1, Arg2, getMinimum(Arg0, Arg1));
    }
    void test_case_2()
    {
        int Arr0[] = {12571295, 2174218, 2015120};
        vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
        int Arg1 = 1;
        int Arg2 = 0;
        verify_case(2, Arg2, getMinimum(Arg0, Arg1));
    }
    void test_case_3()
    {
        int Arr0[] = {5,2,4,52,62,9,8,3,1,11,6};
        vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
        int Arg1 = 11;
        int Arg2 = 3;
        verify_case(3, Arg2, getMinimum(Arg0, Arg1));
    }
    void test_case_4()
    {
        int Arr0[] = {503, 505, 152, 435, 491, 512, 1023, 355, 510, 500, 502, 255, 63, 508, 509, 511, 60, 250, 254, 346};
        vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0])));
        int Arg1 = 510;
        int Arg2 = 5;
        verify_case(4, Arg2, getMinimum(Arg0, Arg1));
    }

// END CUT HERE

};

// BEGIN CUT HERE
int main()
{
    ORSolitaire ___test;
    ___test.run_test(-1);
    return 0;
}
// END CUT HERE
