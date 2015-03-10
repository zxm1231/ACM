/*
     ID: xinming2
     PROG: runround
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
int touth[20];
bool diff(int n)
{
    int a[10] = {0};
    int ok = 1;
    while(n){
        if(a[n % 10] || n % 10 == 0){
            ok = 0;
            break;
        }
        a[n % 10] = 1;
        n /= 10;
    }
    return ok;
}
bool check(int n)
{
    clr(touth , 0);
    string str = "";
    while(n){
        str += n % 10 + '0';
        n /= 10;
    }
    reverse(str.begin() , str.end());
    int l = str.sz;
    str = str + str;
    int pos = str[0];
    int ok = 0;
    for(int i = 0 ; i < 2 * l ; ){
        int num = (i + str[i] - '0') % l;
        if(touth[str[i] - '0']){
            ok = 0;
            break;
        }
        touth[str[i] - '0'] = 1;
        int flag = 1;
        for(int j = 0 ; j < l ; j++){
            if(!touth[str[j] - '0']){
                flag = 0;
                break;
            }
        }
//        cout << i << ' ' << str[i] << ' ' << num << ' ' << flag << endl;
        if(i && str[num] == pos && flag){
            ok = 1;
            break;
        }

        i = num;
    }
    if(ok || l == 1)return 1;
    return 0;
}
int main()
{
    freopen("runround.in" , "r" , stdin);
    freopen("runround.out" , "w" , stdout);
    while(~scanf("%d" , &n)){
//        cout << check(n) << endl;
        int i;
        for(i = n + 1; ; i++){
            if(!diff(i))continue;
//            cout << i  << ' '  << check(i) << endl;
            if(check(i)){
                break;
            }
        }
        printf("%d\n" , i);
    }
    return 0;
}
/*


*/

