/*
     ID: xinming2
     PROG: subset
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

const int MAXN = 60 + 5;
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

LL dp[MAXN][MAXN * MAXN];
int n;
int main()
{
    freopen("subset.in" , "r" , stdin);
    freopen("subset.out" , "w" , stdout);
    while(~scanf("%d" , &n)){
        int sum = n * (n + 1) / 2;
        if(sum & 1){
            puts("0");
        }
        else{
            sum /= 2;
            clr(dp , 0);
            dp[1][0] = dp[1][1] = 1;
            for(int i = 2 ; i <= n ; i++){
                for(int j = 0 ; j <= sum ; j++){
                    dp[i][j] = dp[i - 1][j];
                    if(j >= i)dp[i][j] += dp[i - 1][j - i];
                }
            }
//             for(int i = 0; i <= n; i++) dp[i][0] = 1;
//             for(int i = 1; i <= n; i++){
//                  for(int j = 1; j <= sum; j++){
//                       dp[i][j] += dp[i - 1][j];
//                       if(i <= j) dp[i][j] += dp[i - 1][j - i];
//                  }
//             }
//            for(int i = 1 ; i <= n ; i++){
//                for(int j = 1 ; j <= sum ; j++){
//                    cout << dp[i][j] <<' ' ;
//                }cout << endl;
//            }
            printf("%lld\n" , dp[n][sum] / 2);
        }
    }
    return 0;
}
/*


*/

