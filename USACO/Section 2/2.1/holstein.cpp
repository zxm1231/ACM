/*
     ID: xinming2
     PROG: holstein
     LANG: C++
*/

#include <algorithm>
#include <bitset>
#include <cctype>
#include <cerrno>
#include <clocale>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <limits>
#include <list>
#include <map>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>
#include <cwchar>
#include <cwctype>

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

const int MAXN = 100 + 5;
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

int n , m;
int tar[MAXN];
int has[MAXN];
int a[MAXN][MAXN];
int main()
{
    freopen("holstein.in" , "r" , stdin);
    freopen("holstein.out" , "w" , stdout);
    while(~scanf("%d" , &n)){
        for(int i = 0 ; i < n ; i++){
            scanf("%d" , &tar[i]);
        }
        scanf("%d" , &m);
        for(int i = 0 ; i < m ; i++){
            for(int j = 0 ; j < n ; j++){
                scanf("%d" , &a[i][j]);
            }
        }

        int ans = inf;
        vec ansv;
        for(int pos = 0 ; pos < (1 << m) ; pos++){
            clr(has , 0);
            vec vi;
            int num = 0;
            for(int i = 0 ; i < m ; i++){
                if(pos & (1 << i)){
                    num++;
                    vi.pb(i + 1);
                    for(int j = 0 ; j < n ; j++){
                        has[j] += a[i][j];
                    }
                }
            }
            int ok = 1;
            for(int j = 0 ; j < n ; j++){
                if(has[j] < tar[j]){
                    ok = 0;
                    break;
                }
            }
            if(ok && num < ans){
                ans = num;
                ansv.clear();
                for(int i = 0 ; i < vi.sz ; i++){
                    ansv.pb(vi[i]);
                }
            }
        }
        sort(ansv.begin() , ansv.end());
        printf("%d" , ans);
        for(int i = 0 ; i < ansv.sz ; i++){
            printf(" %d" , ansv[i]);
        }
        puts("");
    }
    return 0;
}
/*


*/


