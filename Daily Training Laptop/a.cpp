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

const int MAXN = 2000 + 5;
const int MAXE = 200000 + 50;
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
map <string , int> mp;
void init()
{
    mp["January"] = 1;
    mp["February"] = 2;
    mp["March"] = 3;
    mp["April"] = 4;
    mp["May"] = 5;
    mp["June"] = 6;
    mp["July"] = 7;
    mp["August"] = 8;
    mp["September"] = 9;
    mp["October"] = 10;
    mp["November"] = 11;
    mp["December"] = 12;

}
bool isLeap(int year)
{
    if(year % 400 == 0)return 1;
    if(year % 4 == 0 && year % 100)return 1;
    return 0;
}
int cmp(int y1 , int m1 , int d1 , int y2 , int m2 , int d2)
{
    if(y1 == y2 && m1 == m2 && d1 == d2)return 2;
    if(y1 == y2){
        if(m1 == m2){
            return d1 < d2;
        }
        return m1 < m2;
    }
    return y1 < y2;
}
int main()
{
    init();
    int _T;
    cin >> _T;
    for(int cas = 1 ; cas <= _T ; cas++){
        int y1 , m1 , d1 , y2,  m2 , d2;
        char str1[20] , str2[20];
        scanf("%s%d,%d" , str1 , &d1 , &y1);
        scanf("%s%d,%d" , str2 , &d2 , &y2);
        m1 = mp[str1];
        m2 = mp[str2];
        int ans = 0;
        if(cmp(y1 , m1 , d1 , y2 , m2 , d2) == 2){
            if(cmp(y1 , m1 , d1 , y1 , 2 , 29) == 2)printf("Case #%d: 1\n" , cas);
            else printf("Case #%d: 0\n" , cas);
        }
        else{
            if(y1 == y2){
                if(isLeap(y1) && cmp(y1 , m1 , d1 , y1 , 2 , 29) && cmp(y2 , 2 , 29 , y2 , m2 , d2))printf("Case #%d: 1\n" , cas);
                else printf("Case #%d: 0\n" , cas);
            }
            else{
                if(isLeap(y1) && cmp(y1 , m1 , d1 , y1 , 2 , 29))ans++;
                if(isLeap(y2) && cmp(y2 , 2 , 29 , y2 , m2 , d2))ans++;
                int a = 0  , b = 0;
                for(int i = y1 + 1 ;; i++){
                    if(isLeap(i)){
                        a = i / 4 - i / 100 + i / 400;
                        break;
                    }
                }
                for(int i = y2 - 1 ;; i--){
                    if(isLeap(i)){
                        b = i / 4 - i / 100 + i / 400;
                        break;
                    }
                }
                ans += b - a + 1;
                printf("Case #%d: %d\n" , cas , ans);
            }

        }
    }
    return 0;
}
/*

*/
