/*
     ID: xinming2
     PROG: sort3
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

const int MAXN = 1000 + 5;
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
int medal[MAXN];
int cntequ(int a , int st , int l)
{
    int num = 0;
    for(int i = st ; i < st + l ; i++){
        if(medal[i] == a)num++;
    }
    return num;
}
int cntnotequ(int a , int st , int l)
{
    int num = 0;
    for(int i = st ; i < st + l ; i++){
        if(medal[i] != a)num++;
    }
    return num;
}
int main()
{
    freopen("sort3.in" , "r" , stdin);
    freopen("sort3.out" , "w" , stdout);
    while(~scanf("%d" , &n)){
        int num1 = 0 , num3 = 0;
        for(int i = 0 ; i < n ; i++){
            scanf("%d" , &medal[i]);
            if(medal[i] == 1)num1++;
            if(medal[i] == 3)num3++;
        }
        int c1 = cntnotequ(1 , 0 , num1);
        int c3 = cntnotequ(3 , n - num3 , num3);

        int c31 = cntequ(3 , 0 , num1);
        int c13 = cntequ(1 , n - num3 , num3);

        int ans = c1 + c3 - min(c13 , c31);
        printf("%d\n" , ans);
    }
    return 0;
}
/*


*/

