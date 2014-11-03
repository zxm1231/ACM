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

const int MAXN = 200 + 5;
const int MAXE = 400000 + 50;
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

int V;
int matchx[MAXN] , matchy[MAXN];
int match;
int vis[MAXN];
int num;
int head[MAXN];
struct Edge
{
    int u , v , next;
}e[MAXE];
void init()
{
    num = 0;
    clr(head , -1);
}
void adde(int u , int v)
{
    e[num].u = u;
    e[num].v = v;
    e[num].next = head[u];
    head[u] = num++;
}
int prehungary()
{
    int res = 0 , u , v;
    for(int i = 1 ; i <= V ; i++){
        u = i;
        for(int j = head[u] ; ~j ; j = e[j].next){
            v = e[j].v;
            if(matchy[v] == -1){
                matchx[u] = v;
                matchy[v] = u;
                res++;
                break;
            }
        }
    }
    return res;
}
int dfs(int u)
{
    int v;
    for(int j = head[u] ; ~j ; j = e[j].next){
        v = e[j].v;
        if(vis[v])continue;
        vis[v] = 1;
        if(matchy[v] == -1 || dfs(matchy[v])){
            matchy[v] = u;
            matchx[u] = v;
            return 1;
        }
    }
    return 0;
}
void solve()
{
    match = 0;
    clr(matchx , -1);
    clr(matchy , -1);
    match += prehungary();
    for(int i = 1 ; i <= V ; i++){
        if(matchx[i] == -1){
            clr(vis , 0);
            if(dfs(i))match++;
        }
    }
}
int main()
{
    int n , m;
    while(~scanf("%d%d" , &n , &m)){
        init();
        int u , v;
        while(scanf("%d%d" , &u , &v) , u > 0){
            adde(u , v);
        }
        V = m;
        solve();
        if(!match){
            puts("No Solution!");
            continue;
        }
        printf("%d\n" , match);
        for(int i = 1 ; i <= n ; i++){
            if(matchx[i] != -1){
                printf("%d %d\n" , i , matchx[i]);
            }
        }
    }
    return 0;
}
