/*
     ID: xinming2
     PROG: hamming
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

int n , b , d;
int cnthammingdist(int x , int y)
{
    bitset <8> ba(x);
    bitset <8> bb(y);
    int num = 0;
    for(int i = 0 ; i < 8 ; i++){
        if(ba[i] != bb[i])num++;
    }
    return num;
}
int main()
{
    freopen("hamming.in" , "r" , stdin);
    freopen("hamming.out" , "w" , stdout);
    while(~scanf("%d%d%d" , &n , &b , &d)){
        vec vi;
        vi.pb(0);
        for(int i = 1 ; i < n ; i++){
            int x = vi[i - 1] + 1;
            while(1){
                int ok = 1;
                for(int j = 0 ; j < i ; j++){
                    if(cnthammingdist(x , vi[j]) < d){
                        ok = 0;
                        break;
                    }
                }
                if(ok)break;
                x++;
            }
            vi.pb(x);
        }
        int i;
        for(i = 0 ; i < n ; i++){
            if(i % 10)putchar(' ');
            printf("%d" , vi[i]);
            if(i % 10 == 9)puts("");
        }
        if(i % 10)puts("");
    }
    return 0;
}
/*


*/


