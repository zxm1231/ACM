/*
     ID: xinming2
     PROG: castle
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

int n , m , cnt;
int U[MAXN][MAXN] , D[MAXN][MAXN] , L[MAXN][MAXN] , R[MAXN][MAXN];
int col[MAXN][MAXN] , SZ[MAXN * MAXN];
void init()
{
    clr(U , 0);clr(D , 0); clr(L , 0);clr(R , 0);
    clr(col , 0); clr(SZ , 0);
}
void flood_fill()
{
    cnt = 0;
    for(int i = 1 ; i <= n ; i++){
        for(int j = 1; j <= m ; j++){
            if(!col[i][j]){
//                cout << "START FROM : " << i << ' ' << j << ' ' << col[i][j] << endl;
                queue <pii> que;
                que.push(mk(i , j));
                col[i][j] = ++cnt;
                SZ[cnt] = 1;
                while(!que.empty()){
//                    outstars
                    pii tmp = que.front();
                    que.pop();
                    int ii = tmp.AA , jj = tmp.BB ;
//                    cout << ii << ' '  << jj << ' ' << cnt << ' ' << endl;
//                    cout << U[ii][jj] << ' ' << D[ii][jj] << ' ' << L[ii][jj] << ' ' << R[ii][jj] << endl;
                    if(col[ii][jj] && (ii != i || jj != j))continue;
                    if(!col[ii][jj]){
                        col[ii][jj] = cnt;
                        SZ[cnt]++;
                    }
                    if(U[ii][jj] && !col[ii - 1][jj])que.push(mk(ii - 1 , jj));
                    if(D[ii][jj] && !col[ii + 1][jj])que.push(mk(ii + 1 , jj));
                    if(L[ii][jj] && !col[ii][jj - 1])que.push(mk(ii , jj - 1));
                    if(R[ii][jj] && !col[ii][jj + 1])que.push(mk(ii , jj + 1));
//                    cout << que.sz << endl;
                }
            }
        }
    }
//    for(int i = 1 ; i <= n ; i++){
//        for(int j = 1 ; j <= m ; j++){
//            cout << col[i][j] << ' ' ;
//        }cout << endl;
//    }
//    for(int i = 1 ; i <= cnt ; i++){
//        cout << i << ' ' << SZ[i] << endl;
//    }
}
int main()
{
    freopen("castle.in" , "r" , stdin);
    freopen("castle.out" , "w" , stdout);
    while(~scanf("%d%d" , &m , &n)){
        init();
        for(int i = 1 ; i <= n ; i++){
            for(int j = 1 ; j <= m ; j++){
                int num;
                scanf("%d" , &num);
                bitset <4> b(num);
                U[i][j] = !b[1];
                D[i][j] = !b[3];
                L[i][j] = !b[0];
                R[i][j] = !b[2];
            }
        }
        flood_fill();

        int maxnum = 0;
        for(int i = 1 ; i <= cnt ; i++){
            maxnum = max(SZ[i] , maxnum);
        }
        printf("%d\n%d\n" , cnt , maxnum);
        ///WSNE
        vector<pair<int , pair <int , pair<int  , int > > > > vp;
        string ansstr = "WSNE";
        for(int i = 1 ; i <= n ; i++){
            for(int j = 1 ; j <= m ; j++){
                if(col[i][j] && col[i - 1][j] && col[i][j] != col[i - 1][j] && !U[i][j]){
//                    cout << i << ' ' << j << "    down   " << endl;
                    vp.pb(mk(n * m - SZ[col[i][j]] - SZ[col[i - 1][j]] , mk(j , mk(n - i , 2))));
                }
                if(col[i][j] && col[i + 1][j] && col[i][j] != col[i + 1][j] && !D[i][j]){
//                    cout << i << ' ' << j << "    down   " << endl;
                    vp.pb(mk(n * m - SZ[col[i][j]] - SZ[col[i + 1][j]] , mk(j , mk(n - i , 1))));
                }
                if(col[i][j] && col[i][j + 1] && col[i][j] != col[i][j + 1] && !R[i][j]){
//                    cout << i << ' ' << j << "   right   " << endl;
                    vp.pb(mk(n * m - SZ[col[i][j]] - SZ[col[i][j + 1]] , mk(j , mk(n - i , 3))));
                }
                if(col[i][j] && col[i][j - 1] && col[i][j] != col[i][j - 1] && !R[i][j]){
//                    cout << i << ' ' << j << "   right   " << endl;
                    vp.pb(mk(n * m - SZ[col[i][j]] - SZ[col[i][j - 1]] , mk(j , mk(n - i , 0))));
                }
            }
        }
        sort(vp.begin() , vp.end());
//        for(int i = 0 ; i < vp.sz ; i++){
//            cout << i << ' ' << n * m - vp[i].AA  << ' ' << n - vp[i].BB.BB.AA << ' ' << vp[i].BB.AA << ' ' << vp[i].BB.BB.BB<< endl;
//        }
//        reverse(vp.begin() , vp.end());
        printf("%d\n%d %d %c\n" , n * m - vp[0].AA , n - vp[0].BB.BB.AA , vp[0].BB.AA , ansstr[vp[0].BB.BB.BB]);
    }
    return 0;
}
/*
32 32
3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 6 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 6
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4
9 8 8 8 8 8 8 8 8 8 8 8 8 8 8 12 9 8 8 8 8 8 8 8 8 8 8 8 8 8 8 12

8 8
3 2 2 6 3 2 2 6
1 0 0 4 1 0 0 4
1 0 0 4 1 0 0 4
1 0 0 4 1 0 0 4
1 0 0 4 1 0 0 4
1 0 0 4 1 0 0 4
1 0 0 4 1 0 0 4
9 8 8 12 9 8 8 12

*/
