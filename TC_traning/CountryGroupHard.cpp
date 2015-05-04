// BEGIN CUT HERE

// END CUT HERE
#line 5 "CountryGroupHard.cpp"
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
#define lson l , mid , rt << 1
#define rson mid + 1 , r , rt << 1 | 1

const int MAXN = 500 + 5;
const int MAXE = 50000 + 50;
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

struct Edge
{
    int v , c , next;
}e[MAXE];
int dis[MAXN];
int in[MAXN];
int head[MAXN];
int cnt[MAXN];
int num , n;
void init()
{
    num = 0;
    clr(head , -1);
}
void adde(int u , int v , int c)
{
//    cout << u << " ---->  " << v << ' ' << c << endl;
    e[num].v = v;
    e[num].c = c;
    e[num].next = head[u];
    head[u] = num++;
}
bool spfa(int s)
{
    queue <int> q;
    clr(cnt , 0);
    clr(in , 0);
    clr(dis , 0x3f);
    dis[s] = 0;
    q.push(s);
    cnt[s] = 1;
    in[s] = 1;
    while(!q.empty()){
        int u = q.front();
        q.pop();
        in[u] = 0;
        for(int j = head[u] ; j != -1 ; j = e[j].next){
            int v = e[j].v;
            if(dis[v] > dis[u] + e[j].c){
                dis[v] = dis[u] + e[j].c;
                if(!in[v]){
                    q.push(v);
                    cnt[v]++;
                    in[v] = 1;
                    if(cnt[v] > n){
                        return 0;
                    }
                }
            }
        }

    }
    return 1;
}
class CountryGroupHard
{
        public:
        string solve(vector <int> a)
        {
            init();
            int ok = 0;
            n = a.sz;
            for(int i = 2  ; i < a.sz - 1; i++){
                if(a[i] && a[i] == a[i - 1]){
                    ok = 1;
                    break;
                }
            }
            if(ok)return "Insufficient";
            int t = 200;
            for(int i = 0 ; i < a.sz ; i++){
                if(!a[i] && i && i != a.sz - 1){
                    if(i == a.sz - 1)adde(i , t , 0);
                    else adde(i , i + 1 , 0);
                }
                else {
                    for(int j = max(i - 1 - a[i] , 0) ; j <= i ; j++){
                        if(j + a[i] > a.sz)continue;
                        if(j + a[i] == a.sz)adde(j , t , 0);
                        else adde(j , j + a[i] , 0);
                    }
                }
            }
            adde(0 , t , inf);
            spfa(0);
//            for(int i = 0 ; i < a.sz ; i++){
//                cout << i << ' ' << dis[i] << endl;
//            }
//            cout << t << ' ' << dis[t] << endl;
            if(!dis[t])return "Sufficient";
            else return "Insufficient";
        }

// BEGIN CUT HERE
	public:
	void run_test(int Case) { if ((Case == -1) || (Case == 0)) test_case_0(); if ((Case == -1) || (Case == 1)) test_case_1(); if ((Case == -1) || (Case == 2)) test_case_2(); if ((Case == -1) || (Case == 3)) test_case_3(); if ((Case == -1) || (Case == 4)) test_case_4(); }
	private:
	template <typename T> string print_array(const vector<T> &V) { ostringstream os; os << "{ "; for (typename vector<T>::const_iterator iter = V.begin(); iter != V.end(); ++iter) os << '\"' << *iter << "\","; os << " }"; return os.str(); }
	void verify_case(int Case, const string &Expected, const string &Received) { cerr << "Test Case #" << Case << "..."; if (Expected == Received) cerr << "PASSED" << endl; else { cerr << "FAILED" << endl; cerr << "\tExpected: \"" << Expected << '\"' << endl; cerr << "\tReceived: \"" << Received << '\"' << endl; } }
	void test_case_0() { int Arr0[] = {0,2,3,0,0}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); string Arg1 = "Sufficient"; verify_case(0, Arg1, solve(Arg0)); }
	void test_case_1() { int Arr0[] = {0,2,0}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); string Arg1 = "Insufficient"; verify_case(1, Arg1, solve(Arg0)); }
	void test_case_2() { int Arr0[] = {0,3,0,0,3,0}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); string Arg1 = "Sufficient"; verify_case(2, Arg1, solve(Arg0)); }
	void test_case_3() { int Arr0[] = {0,0,3,3,0,0}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); string Arg1 = "Insufficient"; verify_case(3, Arg1, solve(Arg0)); }
	void test_case_4() { int Arr0[] = {2,2,0,2,2}; vector <int> Arg0(Arr0, Arr0 + (sizeof(Arr0) / sizeof(Arr0[0]))); string Arg1 = "Sufficient"; verify_case(4, Arg1, solve(Arg0)); }

// END CUT HERE

};

// BEGIN CUT HERE
int main()
{
        CountryGroupHard ___test;
        ___test.run_test(-1);
        return 0;
}
// END CUT HERE
