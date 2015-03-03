/*
     ID: xinming2
     PROG: wormhole
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


int n;

int X[MAXN] , Y[MAXN];
int right_p[MAXN];
int par[MAXN];

bool circle_exist()
{
    for(int i = 1 ; i <= n ; i++){
        int pos = i;
        for(int cnt = 1 ; cnt <= n ; cnt++){
            pos = right_p[par[pos]];
        }
        if(pos)return 1;
    }
    return 0;
}
int solve()
{
    int i , tot = 0;
    for(i = 1 ; i <= n ; i++){
        if(!par[i])break;
    }
//    cout << circle_exist() << endl;
//    cout << i << endl;
    if(i > n)return circle_exist();
    for(int j = i + 1 ; j <= n ; j++){
        if(!par[j]){
            par[j] = i;
            par[i] = j;
            tot += solve();
            par[i] = par[j] = 0;
        }
    }
    return tot;
}
int main()
{
    freopen("wormhole.in" , "r" , stdin);
    freopen("wormhole.out" , "w" , stdout);
    while(~scanf("%d" , &n)){
        for(int i = 1 ; i <= n ; i++){
            scanf("%d%d" , &X[i] , &Y[i]);
        }
        clr(right_p , 0);
        clr(par , 0);
        for(int i = 1 ; i <= n ; i++){
            for(int j = 1 ; j <= n ; j++){
                if(X[i] < X[j] && Y[i] == Y[j]){
                    if(!right_p[i] || (X[j] - X[i]) < (X[right_p[i]] - X[i])){
                        right_p[i] = j;
                    }
                }
            }
        }
//        for(int i = 1 ; i <= n ; i++){
//            cout << X[i] << ' ' << Y[i] << " _____ " << X[right_p[i]] << ' '<< Y[right_p[i]] << endl;
//            cout << i << ' ' << right_p[i] <<endl;
//        }
        printf("%d\n" , solve());
    }
    return 0;
}
/*


*/
