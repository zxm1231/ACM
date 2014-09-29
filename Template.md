hu#配置环境

- codeblocks    在codeblocks --> setting -->环境变量

		terminal： gnome-terminal -t $TITLE -x   
- eclipse

	Eclipse中默认是输入"."后出现自动提示，用于类成员的自动提示，可是有时候我们希望它能在我们输入类的首字母后就出现自动提示，可以节省大量的输入时间（虽然按alt + /会出现提示，但还是要多按一次按键，太麻烦了）。    从Window -> preferences -> Java -> Editor -> Content assist -> Auto-Activation下，我们可以在"."号后面加入我们需要自动提示的首字幕，比如"ahiz"。    然后我们回到Eclipse的开发环境，输入"a"，提示就出现了。但是我们可以发现，这个Auto-Activation下的输入框里最多只能输入5个字母，也许是Eclipse的开发人员担心我们输入的太多会影响性能，但计算机的性能不用白不用，所以我们要打破这个限制。其实上面都是铺垫，制造一下气氛，以显得我们下面要做的事情很牛似的，其实不然，一切都很简单。嘿嘿 :)
	在"."后面随便输入几个字符，比如"abij"，然后回到开发环境，File -> export -> general -> preferences -> 选一个地方保存你的首选项，比如C:"a.epf用任何文本编辑器打开a.epf，查找字符串“abij”，找到以后，替换成“abcdefghijklmnopqrstuvwxyz”，总之就是你想怎样就怎样！！然后回到Eclipse，File -> import -> general -> preferences -> 导入刚才的a.epf文件。此时你会发现输入任何字幕都可以得到自动提示了。爽！！！最后：自动提示弹出的时间最好改成100毫秒以下，这样会比较爽一点，不然你都完事了，自动提示才弹出来:)，不过也要看机器性能。
	FileWriter fileWriter=new FileWriter("c:\\Result.txt");
	int [] a=new int[]{11112,222,333,444,555,666};
	for (int i = 0; i < a.length; i++) {
	fileWriter.write(String.valueOf(a[i])+" ");
	}
##头文件

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
	//#define Online_Judge
	#define outstars cout << "***********************" << endl;
	#define clr(a,b) memset(a,b,sizeof(a))
	#define lson l , mid  , rt << 1
	#define rson mid + 1 , r , rt << 1 | 1
	#define mk make_pair
	#define pb push_back
	#define sz size()
	#define AA first
	#define BB second
	#define eps 1e-10
	#define zero(x) fabs(x) < eps
	#define deq(a , b) fabs(a - b) < eps
	
	const int MAXN = 40000 + 50;
	const int MAXV = 4000 + 50;
	const int MAXQ = 1000000 + 50;
	const int sigma_size = 26;
	 
	const int inf = 0x3f3f3f3f;
	const int INF = ~0U >> 1;
	const long long LLinf = 0x3FFFFFFFFFFFFFFFLL;
	const long long LLINF = 1 << 63 - 1;
	const int IMIN = 0x80000000;
	 
	 
	const long long MOD = 1000000000 + 7;
	const int mod = 10007;
	typedef long long LL;
	const double PI = acos(-1.0);
	typedef pair<int , int> pii;
	typedef vector<int> vec;
	typedef vector<vec> mat;
	 
	#define Bug(s) cout << "s = " << s << endl;
	///#pragma comment(linker, "/STACK:102400000,102400000")
	 
	int main()
	{
	    //ios::sync_with_stdio(false);
	    #ifdef Online_Judge
	        freopen(".in","r",stdin);
	        freopen(".out","w",stdout);
	    #endif Online_Judge
	    
	    return 0;
	}
##输入输出挂
	int Scan()
	{
	    int flag = 1;
	    char ch;
	    int a = 0;
	    while((ch = getchar()) == ' ' || ch == '\n');
	    if(ch == '-') flag = -1;
	    else
	    a += ch - '0';
	    while((ch = getchar()) != ' ' && ch != '\n')
	    {
	        a *= 10;
	        a += ch - '0';
	    }
	    return flag * a;
	}
	void Out(int a)
	{
	    if(a < 0) {putchar('-'); a = -a;}
	    if(a >= 10)Out(a / 10);
	    putchar(a % 10 + '0');
	}

#杂项
##C++大整数类

	#define MAXN 9999
    #define MAXSIZE 10
    #define DLEN 4 
	class BigNum
	{ 
	private: 
		int a[500];    //可以控制大数的位数 
		int len;       //大数长度
	public: 
		BigNum(){ len = 1;memset(a,0,sizeof(a)); }   //构造函数
		BigNum(const int);       //将一个int类型的变量转化为大数
		BigNum(const char*);     //将一个字符串类型的变量转化为大数
		BigNum(const BigNum &);  //拷贝构造函数
		BigNum &operator=(const BigNum &);   //重载赋值运算符，大数之间进行赋值运算
	 
		friend istream& operator>>(istream&,  BigNum&);   //重载输入运算符
		friend ostream& operator<<(ostream&,  BigNum&);   //重载输出运算符
	 
		BigNum operator+(const BigNum &) const;   //重载加法运算符，两个大数之间的相加运算 
		BigNum operator-(const BigNum &) const;   //重载减法运算符，两个大数之间的相减运算 
		BigNum operator*(const BigNum &) const;   //重载乘法运算符，两个大数之间的相乘运算 
		BigNum operator/(const int   &) const;    //重载除法运算符，大数对一个整数进行相除运算
	 
		BigNum operator^(const int  &) const;    //大数的n次方运算
		int    operator%(const int  &) const;    //大数对一个int类型的变量进行取模运算    
		bool   operator>(const BigNum & T)const;   //大数和另一个大数的大小比较
		bool   operator>(const int & t)const;      //大数和一个int类型的变量的大小比较
	 
		void print();       //输出大数
	}; 
	BigNum::BigNum(const int b)     //将一个int类型的变量转化为大数
	{ 
		int c,d = b;
		len = 0;
		memset(a,0,sizeof(a));
		while(d > MAXN)
		{
			c = d - (d / (MAXN + 1)) * (MAXN + 1); 
			d = d / (MAXN + 1);
			a[len++] = c;
		}
		a[len++] = d;
	}
	BigNum::BigNum(const char*s)     //将一个字符串类型的变量转化为大数
	{
		int t,k,index,l,i;
		memset(a,0,sizeof(a));
		l=strlen(s);   
		len=l/DLEN;
		if(l%DLEN)
			len++;
		index=0;
		for(i=l-1;i>=0;i-=DLEN)
		{
			t=0;
			k=i-DLEN+1;
			if(k<0)
				k=0;
			for(int j=k;j<=i;j++)
				t=t*10+s[j]-'0';
			a[index++]=t;
		}
	}
	BigNum::BigNum(const BigNum & T) : len(T.len)  //拷贝构造函数
	{ 
		int i; 
		memset(a,0,sizeof(a)); 
		for(i = 0 ; i < len ; i++)
			a[i] = T.a[i]; 
	} 
	BigNum & BigNum::operator=(const BigNum & n)   //重载赋值运算符，大数之间进行赋值运算
	{
		int i;
		len = n.len;
		memset(a,0,sizeof(a)); 
		for(i = 0 ; i < len ; i++) 
			a[i] = n.a[i]; 
		return *this; 
	}
	istream& operator>>(istream & in,  BigNum & b)   //重载输入运算符
	{
		char ch[MAXSIZE*4];
		int i = -1;
		in>>ch;
		int l=strlen(ch);
		int count=0,sum=0;
		for(i=l-1;i>=0;)
		{
			sum = 0;
			int t=1;
			for(int j=0;j<4&&i>=0;j++,i--,t*=10)
			{
				sum+=(ch[i]-'0')*t;
			}
			b.a[count]=sum;
			count++;
		}
		b.len =count++;
		return in;
	 
	}
	ostream& operator<<(ostream& out,  BigNum& b)   //重载输出运算符
	{
		int i;  
		cout << b.a[b.len - 1]; 
		for(i = b.len - 2 ; i >= 0 ; i--)
		{ 
			cout.width(DLEN); 
			cout.fill('0'); 
			cout << b.a[i]; 
		} 
		return out;
	}
	 
	BigNum BigNum::operator+(const BigNum & T) const   //两个大数之间的相加运算
	{
		BigNum t(*this);
		int i,big;      //位数   
		big = T.len > len ? T.len : len; 
		for(i = 0 ; i < big ; i++) 
		{ 
			t.a[i] +=T.a[i]; 
			if(t.a[i] > MAXN) 
			{ 
				t.a[i + 1]++; 
				t.a[i] -=MAXN+1; 
			} 
		} 
		if(t.a[big] != 0)
			t.len = big + 1; 
		else
			t.len = big;   
		return t;
	}
	BigNum BigNum::operator-(const BigNum & T) const   //两个大数之间的相减运算 
	{  
		int i,j,big;
		bool flag;
		BigNum t1,t2;
		if(*this>T)
		{
			t1=*this;
			t2=T;
			flag=0;
		}
		else
		{
			t1=T;
			t2=*this;
			flag=1;
		}
		big=t1.len;
		for(i = 0 ; i < big ; i++)
		{
			if(t1.a[i] < t2.a[i])
			{ 
				j = i + 1; 
				while(t1.a[j] == 0)
					j++; 
				t1.a[j--]--; 
				while(j > i)
					t1.a[j--] += MAXN;
				t1.a[i] += MAXN + 1 - t2.a[i]; 
			} 
			else
				t1.a[i] -= t2.a[i];
		}
		t1.len = big;
		while(t1.a[len - 1] == 0 && t1.len > 1)
		{
			t1.len--; 
			big--;
		}
		if(flag)
			t1.a[big-1]=0-t1.a[big-1];
		return t1; 
	} 
	 
	BigNum BigNum::operator*(const BigNum & T) const   //两个大数之间的相乘运算 
	{ 
		BigNum ret; 
		int i,j,up; 
		int temp,temp1;   
		for(i = 0 ; i < len ; i++)
		{ 
			up = 0; 
			for(j = 0 ; j < T.len ; j++)
			{ 
				temp = a[i] * T.a[j] + ret.a[i + j] + up; 
				if(temp > MAXN)
				{ 
					temp1 = temp - temp / (MAXN + 1) * (MAXN + 1); 
					up = temp / (MAXN + 1); 
					ret.a[i + j] = temp1; 
				} 
				else
				{ 
					up = 0; 
					ret.a[i + j] = temp; 
				} 
			} 
			if(up != 0) 
				ret.a[i + j] = up; 
		} 
		ret.len = i + j; 
		while(ret.a[ret.len - 1] == 0 && ret.len > 1)
			ret.len--; 
		return ret; 
	} 
	BigNum BigNum::operator/(const int & b) const   //大数对一个整数进行相除运算
	{ 
		BigNum ret; 
		int i,down = 0;   
		for(i = len - 1 ; i >= 0 ; i--)
		{ 
			ret.a[i] = (a[i] + down * (MAXN + 1)) / b; 
			down = a[i] + down * (MAXN + 1) - ret.a[i] * b; 
		} 
		ret.len = len; 
		while(ret.a[ret.len - 1] == 0 && ret.len > 1)
			ret.len--; 
		return ret; 
	}
	int BigNum::operator %(const int & b) const    //大数对一个int类型的变量进行取模运算    
	{
		int i,d=0;
		for (i = len-1; i>=0; i--)
		{
			d = ((d * (MAXN+1))% b + a[i])% b;  
		}
		return d;
	}
	BigNum BigNum::operator^(const int & n) const    //大数的n次方运算
	{
		BigNum t,ret(1);
		int i;
		if(n<0)
			exit(-1);
		if(n==0)
			return 1;
		if(n==1)
			return *this;
		int m=n;
		while(m>1)
		{
			t=*this;
			for( i=1;i<<1<=m;i<<=1)
			{
				t=t*t;
			}
			m-=i;
			ret=ret*t;
			if(m==1)
				ret=ret*(*this);
		}
		return ret;
	}
	bool BigNum::operator>(const BigNum & T) const   //大数和另一个大数的大小比较
	{ 
		int ln; 
		if(len > T.len)
			return true; 
		else if(len == T.len)
		{ 
			ln = len - 1; 
			while(a[ln] == T.a[ln] && ln >= 0)
				ln--; 
			if(ln >= 0 && a[ln] > T.a[ln])
				return true; 
			else
				return false; 
		} 
		else
			return false; 
	}
	bool BigNum::operator >(const int & t) const    //大数和一个int类型的变量的大小比较
	{
		BigNum b(t);
		return *this>b;
	}
	 
	void BigNum::print()    //输出大数
	{ 
		int i;   
		cout << a[len - 1]; 
		for(i = len - 2 ; i >= 0 ; i--)
		{ 
			cout.width(DLEN); 
			cout.fill('0'); 
			cout << a[i]; 
		} 
		cout << endl;
	}
	int main(void)
	{
		int i,n;
		BigNum x[101];      //定义大数的对象数组
		x[0]=1;
		for(i=1;i<101;i++)
			x[i]=x[i-1]*(4*i-2)/(i+1);
		while(scanf("%d",&n)==1 && n!=-1)
		{
			x[n].print();
		}
	}
##矩阵快速幂
	mat mul(mat & A , mat & B)  
	{  
	    mat C(A.size() , vec(B[0].size()));  
	    for(int i = 0 ; i < A.size() ; i++)  
	    {  
	        for(int k = 0 ; k < B.size() ; k++)  
	        {  
	            for(int j = 0  ; j < B[0].size() ; j++)  
	            {  
	                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % m;  
	            }  
	        }  
	    }  
	    return C;  
	}  
	mat pow(mat A , LL n)  
	{  
	    mat B(A.size() , vec(A.size()));  
	    for(int i = 0 ; i < A.size() ; i++)B[i][i] = 1;  
	    while(n > 0)  
	    {  
	        if(n & 1)B = mul(B , A);  
	        A = mul(A , A);  
	        n >>= 1;  
	    }  
	    return B;  
	} 
##RMQ

###一维

	int a[MAXN];
	int mm[MAXN];
	int dp[MAXN][20];
	int n , k;
	void initRMQ()
	{
		mm[0] = -1;
	    for(int i = 1 ; i <= MAXN ; i++){
	        mm[i] = ((i & (i - 1)) == 0) ? mm[i - 1] + 1 : mm[i - 1];
	    }
	    for(int i = 1 ; i <= n ; i++){
	        dp[i][0] = a[i];
	    }
	    for(int i = 1 ; i <= mm[n] + 1 ; i++){
	        for(int j = 1 ; j + (1 << i) - 1 <= n ; j++){
	            dp[j][i] = max(dp[j][i - 1] , dp[j + (1 << (i - 1))][i - 1]);
	        }
	    }
	}
	int rmq(int L , int R)
	{
	    int k = mm[R - L + 1];
	    return max(dp[L][k] , dp[R - (1 << k) + 1][k]);
	}

###二维
	int a[MAXN][MAXN];
	int mm[MAXN];
	int dp[MAXN][MAXN][9][9];
	int n , m , q;
	void initRMQ()
	{
		mm[0] = -1;
	    for(int i = 1 ; i <= MAXN ; i++){
	        mm[i] = ((i & (i - 1)) == 0) ? mm[i - 1] + 1 : mm[i - 1];
	    }
	    for(int i = 1 ; i <= n ; i++){
	        for(int j = 1 ; j <= m ; j++){
	            dp[i][j][0][0] = a[i][j];
	        }
	    }
	    for(int ii = 0 ; ii <= mm[n] ; ii++){
	        for(int jj = 0 ; jj <= mm[m] ; jj++){
	            if(ii + jj){
	                for(int i = 1 ; i + (1 << ii) - 1 <= n ; i++){
	                    for(int j = 1 ; j + (1 << jj) - 1 <= m ; j++){
	                        if(ii)dp[i][j][ii][jj] = max(dp[i][j][ii - 1][jj] , dp[i + (1 << (ii - 1))][j][ii - 1][jj]);
	                        else dp[i][j][ii][jj] = max(dp[i][j][ii][jj - 1] , dp[i][j + (1 << (jj - 1))][ii][jj - 1]);
	                    }
	                }
	            }
	        }
	    }
	}
	int rmq(int x1 , int y1 , int x2 , int y2)
	{
	    int k1 = mm[x2 - x1 + 1];
	    int k2 = mm[y2 - y1 + 1];
	    x2 = x2 - (1 << k1) + 1;
	    y2 = y2 - (1 << k2) + 1;
	    return max(max(dp[x1][y1][k1][k2] , dp[x1][y2][k1][k2]) , max(dp[x2][y1][k1][k2] , dp[x2][y2][k1][k2]));
	}
##LCA
###离线——Tarjan
	//flag初始化为0 ， 需要并查集的find ， g记录每个点的query对点 ， dis记录该点距根的距离
	//初始tarjan根（1）
	void tarjan(int u)
	{
	    flag[u] = 1;
	    for(int i = 0 ; i < g[u].sz ; i++){
	        int v = g[u][i].v , id = g[u][i].id;
	        if(flag[v])query[id][2] = find(v);
	    }
	    for(int i = head[u] ; i != -1 ; i = e[i].next){
	        int v = e[i].v , w = e[i].c;
	        if(!flag[v]){
	            dis[v] = dis[u] + w;
	            tarjan(v);
	            fa[v] = u;
	        }
	    }
	}

###在线——基于RMQ
	vector<pair<int,int> >e[MAXN];
	int dist[MAXN];
	int root;
	int depth,b[MAXN*2],a[MAXN*2],tot;
	int p[MAXN],f[MAXN];
	int dp[MAXN*2][20];
	int ori[MAXN];
	void init()
	{
	    tot = 0;
	    depth = 0;
	    clr(ori , 0);
	    for(int i = 0 ; i <= n ; i++)e[i].clear();
	}
	int find_root()
	{
	    for(int i = 1 ; i <= n; i++){
	        if(!ori[i])return i;
	    }
	}
	void dfs(){
	    stack<pair<int,pair<int,int> > >s;
	    s.push(mk(root,mk(0,0)));
	    while(!s.empty()){
	        pair<int,pair<int,int> >now=s.top();s.pop();
	        int u=now.first,pre=now.second.first,i=now.second.second;
	        if(i==0){
	            int t=++depth;
	            b[++tot]=t;
	            f[t]=u;
	            p[u]=tot;
	        }
	        if(i<e[u].size()){
	            int v=e[u][i].first,w=e[u][i].second;
	            s.push(mk(u,mk(pre,i+1)));
	            if(v==pre) continue;
	            dist[v]=dist[u]+w;
	            s.push(mk(v,mk(u,0)));
	        }
	        else
	            b[++tot]=b[p[pre]];
	    }
	}
	
	void Init_rmq(int n){
	    for(int i=1;i<=n;i++)
	        dp[i][0]=b[i];
	    int m=floor(log(n*1.0)/log(2.0));
	    for(int j=1;j<=m;j++)
	        for(int i=1;i<=n-(1<<j)+1;i++)
	            dp[i][j]=min(dp[i][j-1],dp[i+(1<<(j-1))][j-1]);
	}
	int rmq(int l,int r){
	    int k=floor(log((r-l+1)*1.0)/log(2.0));
	    return min(dp[l][k],dp[r-(1<<k)+1][k]);
	}
	int lca(int a,int b){
	    if(p[a]>p[b]) swap(a,b);
	    return f[rmq(p[a],p[b])];
	}

##最小覆盖圆

#图论
##最短路
###dij
	int dij(){
	    priority_queue<pii, vector<pii>, greater<pii> >Q;
	    memset(vis, false, sizeof(vis));
	    for(int i = 0; i <= 4 * n * m + 2; i++) dis[i] = inf;
	    dis[st] = 0;
	    Q.push(pii(dis[st], st));
	    while(!Q.empty())
	    {
	        pii x = Q.top();
	        Q.pop();
	        int u = x.second;
	        if(vis[u]) continue;
	        vis[u] = true;
	        for(int i = head[u]; i != -1; i = e[i].next)
	        {
	            int v = e[i].v, w = e[i].c;
	            if(dis[v] > dis[u] + w)
	            {
	                dis[v] = dis[u] + w;
	                Q.push(pii(dis[v], v));
	            }
	        }
	    }
	    return dis[en];
	}
###最小环
	///inf不要开太大
	void init()
	{
	    for(int i = 1; i <= n ; i++){
	        for(int j = 1 ; j <= n ;j ++){
	            g[i][j] = dist[i][j] = inf;
	            pre[i][j] = i;
	        }
	    }
	}
	void floyd()
	{
	    for(int k = 1 ; k <= n ; k++){
			//求最大环不要这部分循环, 普通的floyd就可以
			//最后最大环为max(dist[i][i])
	        for(int i = 1;  i < k ; i++){
	            for(int j = i + 1 ; j < k ; j++){
	                if(ans > g[i][j] + dist[i][k] + dist[k][j]){
	                    ans = g[i][j] + dist[i][k] + dist[k][j];
	                    num = 0;
	                    int p = j;
	                    while(p != i){
	                        path[num++] = p;
	                        p = pre[i][p];
	                    }
	                    path[num++] = i;
	                    path[num++] = k;
	                }
	            }
	        }
	        for(int i = 1 ; i <= n ; i++){
	            for(int j = 1 ; j <= n ; j++){
	                if(g[i][j] > g[i][k] + g[k][j]){
	                    g[i][j] = g[i][k] + g[k][j];
	                    pre[i][j] = pre[k][j];
	                }
	            }
	        }
	    }
	}

###最短路次短路计数
	void dij()
	{
	    clr(cnt , 0);
	    clr(vis , 0);
	    for(int i = 0 ; i <= n ; i++){
	        dis[i][0] = inf;
	        dis[i][1] = inf;
	    }
	    dis[st][0] = 0;
	    cnt[st][0] = 1;
	    for(int i = 1 ; i <= 2 * n ; i++){
	        int tmp = inf , tx = -1 , flag;
	        for(int j = 1 ; j <= n ; j++){
	            if(!vis[j][0] && dis[j][0] < tmp){
	                tmp = dis[j][0];
	                tx = j;
	                flag = 0;
	            }
	            else if(!vis[j][1] && dis[j][1] < tmp){
	                tmp = dis[j][1];
	                tx = j;
	                flag = 1;
	            }
	        }
	        if(tx == -1)break;
	        vis[tx][flag] = 1;
	        for(int j = head[tx] ; ~j ; j = e[j].next){
	            int v = e[j].v , c = e[j].c;
	            if(dis[tx][flag] + c < dis[v][0]){
	                dis[v][1] = dis[v][0];
	                cnt[v][1] = cnt[v][0];
	                dis[v][0] = dis[tx][flag] + c;
	                cnt[v][0] = cnt[tx][flag];
	            }
	            else if(dis[tx][flag] + c == dis[v][0]){
	                cnt[v][0] += cnt[tx][flag];
	            }
	            else if(dis[tx][flag] + c == dis[v][1]){
	                cnt[v][1] += cnt[tx][flag];
	            }
	            else if(dis[tx][flag] + c < dis[v][1]){
	                dis[v][1] = dis[tx][flag] + c;
	                cnt[v][1] = cnt[tx][flag];
	            }
	        }
	    }
	}
##最小生成树
###最小树形图

	最小树形图（根固定） O(VE)
	有向图最小生成树
	根不固定，添加一个根节点与所有点连无穷大的边！
	根据pre的信息能构造出这棵树！
	int n , m;
	int pre[MAXN] , in[MAXN] , vis[MAXN];
	int id[MAXN];
	int Directed_mst(int root)
	{
	    n++;
	    int ret = 0;
	    while(1){
	        for(int i = 0 ; i < n ; i++){
	            in[i]  = inf;
	            id[i] = -1;
	            vis[i] = -1;
	        }
	        for(int i = 0 ; i < m ; i++){
	            int u = e[i].u;
	            int v = e[i].v;
	            if(e[i].c >= in[v] || u == v)continue;
	            pre[v] = u;
	            in[v] = e[i].c;
	        }
	        in[root] = 0;
	        pre[root] = root;
	        for(int i = 0 ; i < n ; i++){
	            ret += in[i];
	            if(in[i] == inf){
	                return -1;
	            }
	        }
	        int cnt = 0;
	        for(int i = 0 ; i < n ; i++){
	            if(vis[i] == -1){
	                int v = i;
	                while(vis[v] == -1){
	                    vis[v] = i;
	                    v = pre[v];
	                }
	                if(vis[v] != i || v  == root)continue;
	                for(int u = pre[v] ; u != v ; u = pre[u]){
	                    id[u] = cnt;
	                }
	                id[v] = cnt++;
	            }
	        }
	        if(!cnt)break;
	        for(int i = 0 ; i < n ; i++){
	            if(id[i] == -1)id[i] = cnt++;
	        }
	        for(int i = 0 ; i < m ; i++){
	            int v = e[i].v;
	            e[i].u = id[e[i].u];
	            e[i].v = id[e[i].v];
	            e[i].c -= in[v];
	        }
	        n = cnt;
	        root = id[root];
	    }
	    return ret;
	}
###斯坦纳树
###度限制最小生成树
###最优比例生成树
	double prim(double s)
	{
	    clr(vis , 0);
	    double retdis = 0;
	    int retcost = 0;
	    vis[1] =  1;
	    for(int i = 2 ; i <= n ; i++){
	        lowcost[i] = abs(v[1].h - v[i].h) - s * dis[1][i];
	        closest[i] = 1;
	    }
	    for(int num = 0 ; num < n - 1 ; num++){
	        double minans = inf;
	        int u;
	        for(int i = 1 ; i <= n ; i++){
	            if(vis[i])continue;
	            if(lowcost[i] < minans){
	                minans = lowcost[i];
	                u = i;
	            }
	        }
	        vis[u] = 1;
	        retcost += abs(v[closest[u]].h - v[u].h);
	        retdis += dis[closest[u]][u];
	        for(int i = 1 ; i <= n ; i++){
	            double cal = abs(v[u].h - v[i].h) - s * dis[u][i];
	            if(!vis[i] &&cal < lowcost[i]){
	                lowcost[i] = cal;
	                closest[i] = u;
	            }
	        }
	    }
	    return 1.0 * retcost / retdis;
	}
	int main()
	{
	    double ans = 0;
        while(1){
        	double tmp = prim(ans);
            if(fabs(tmp - ans) < eps){
                break;
            }
            ans = tmp;
        }
        printf("%.3f\n" , ans);
	}

###最小生成树计数
###生成树计数
##连通性
###SCC	
	注释部分为判断仙人掌图
	void tarjan(int u)
	{
	    low[u] = dfn[u] = ++dindex;
	    st[++tail] = u;
	    instack[u] = 1;
	//  int cnt = 0;
	    for(int j = head[u] ; ~j ; j = e[j].next) {
	        int v = e[j].v;
	//        if(color[v])ok = 0;
	        if(!dfn[v]) {
	            tataarrjan(v);
	//            if(low[v] > dfn[u])ok = 0;
	//            if(low[v] < dfn[u])cnt++;
	//            if(cnt == 2)ok = 0;
	            low[u] = min(low[u] , low[v]);
	        }
	        else if(instack[v]) {
	            low[u] = min(low[u] , dfn[v]);
	//            cnt++;
	//            if(cnt == 2)ok = 0;
	        }
	    }
	    if(dfn[u] == low[u]) {
	        int i;
	        bcnt++;
	        do {
	            i = st[tail--];
	            instack[i] = 0;
	            cmp[i] = bcnt;
	        } while(i != u);
	    }
	//    color[u] = 1;
	}
	void scc()
	{
	    clr(dfn , 0);
	    clr(low , 0);
	    clr(instack , 0);
	    clr(cmp , 0);
	//  clr(color , 0);
	    tail = bcnt = dindex = 0;
	    for(int i = 1 ; i <= N ; i++) {
	        if(!dfn[i])tarjan(i);
	    }
	}

###点联通
###边双联通分量

	void tarjan(int u, int pre)
	{
	    dfn[u] = low[u] = nindex++;
	    instack[u] = 1;
	    st[++top] = u;
	    int v;
	    for(int j = head[u]; ~j ; j = e[j].next) {
	        v = e[j].v;
	        if(v == pre)continue;
	        if(!dfn[v]) {
	            tarjan(v, u);
	            low[u] = min(low[u] , low[v]);
	            if(low[v] > dfn[u]){
	                bridge++;
	                e[j].cut = 1;
	                e[j ^ 1].cut = 1;
	            }
	        }
	        else if(instack[v]) {
	            low[u] = min(low[u] , dfn[v]);
	        }
	    }
	    if(dfn[u] == low[u]) {
	        ncnt++;
	        do {
	            v = st[top--];
	            instack[v] = 0;
	            cmp[v] = ncnt;
	        } while(v != u);
	    }
	}
	void solve()
	{
	    clr(dfn , 0);
	    clr(low , 0);
	    clr(instack , 0);
	    clr(cmp , 0);
	    ncnt = nindex = top = 0;
	    bridge = 0;
	    tarjan(1 , -1);
	}

##匹配
###最大匹配 匈牙利
	int V;
	int match;
	int matchx[MAXN], matchy[MAXN];
	int vis[MAXN];
	int preHungary()
	{
		int res = 0, u, v;
		for(int i = 1; i <= V; i++)
		{
			u = i;
			for(int j = head1[i]; ~j ; j = e1[j].next)
			{
				v = e1[j].v;
				if(matchy[v] == -1)
				{
					matchx[u] = v;
					matchy[v] = u;
					res++;
					break;
				}
			}
		}
		return res;
	}
	bool dfs(int u)
	{
		int v;
		for(int j = head1[u]; ~j ; j = e1[j].next)
		{
			v = e1[j].v;
			if(vis[v]) continue;
			vis[v] = 1;
			if(matchy[v] == -1 || dfs(matchy[v]))
			{
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
		clr(matchx, -1); clr(matchy, -1);
		match += preHungary();
		for(int i = 1; i <= V; i++)
			if(matchx[i] == -1)
			{
				clr(vis, 0);
				if(dfs(i)) match++;
			}
	}

##网络流
###最大流
####Dinic
	void adde(int u , int v , int c , int cc)
	{
	    e[num].u = u , e[num].v = v , e[num].c = c , e[num].next = head[u] , head[u] = num++;
	    e[num].u = v , e[num].v = u , e[num].c = cc , e[num].next = head[v] , head[v] = num++;
	}
	int dis[MAXN] , cur[MAXN] , sta[MAXN] , que[MAXN] , pre[MAXN];
	bool bfs(int s , int t , int n)
	{
	    int front = 0 , tail = 0;
	    clr(dis , -1);
	    dis[s] = 0;
	    que[tail++] = s;
	    while(front < tail){
	        for(int i = head[que[front ++ ]] ; i != -1 ; i = e[i].next){
	            if(e[i].c > 0 && dis[e[i].v] == -1){
	                dis[e[i].v] = dis[e[i].u] + 1;
	                if(e[i].v == t)return 1;
	                que[tail ++] = e[i].v;
	            }
	        }
	    }
	    return 0;
	}
	int dinic(int s , int t , int n)
	{
	    int maxflow = 0;
	    while(bfs(s , t , n)){
	        for(int i = 0 ; i < n ; i++)cur[i] = head[i];
	        int u = s , tail = 0;
	        while(cur[s] != -1){
	            if(u == t){
	                int det = INF;
	                for(int i = tail - 1 ; i >= 0 ; i--){
	                    det = min(det , e[sta[i]].c);
	                }
	                maxflow += det;
	                for(int i = tail - 1 ; i >= 0 ; i--){
	                    e[sta[i]].c -= det;
	                    e[sta[i] ^ 1].c += det;
	                    if(e[sta[i]].c == 0)tail = i;
	                }
	                u = e[sta[tail]].u;
	            }
	            else if(cur[u] != -1 && e[cur[u]].c > 0 && dis[u] + 1 == dis[e[cur[u]].v]){
	                sta[tail ++] = cur[u];
	                u = e[cur[u]].v;
	            }
	            else{
	                while(u != s && cur[u] == -1){
	                    u = e[sta[--tail]].u;
	                }
	                cur[u] = e[cur[u]].next;
	            }
	        }
	    }
	    return maxflow;
	}
###ISAP
	void adde(int u , int v , int c , int cc)
	{
	    e[num].u = u , e[num].v = v , e[num].c = c , e[num].next = head[u] , head[u] = num++;
	    e[num].u = v , e[num].v = u , e[num].c = cc , e[num].next = head[v] , head[v] = num++;
	}
	 
	int dis[MAXN] , pre[MAXN] , cur[MAXN] , gap[MAXN];
	int ISAP(int s , int t , int n)
	{
	    clr(dis , 0);
	    clr(gap , 0);
	    clr(cur , 0);
	    for(int i = 0 ; i < n ; i++)cur[i] = head[i];
	    gap[0] = n;
	    pre[s] = s;
	    int u = s , maxflow = 0;
	    while(dis[s] <= n){
	        bool flag = false;
	        for(int i = cur[u] ; i != -1;  i = e[i].next){
	            if(e[i].c > 0 && dis[u] == dis[e[i].v] + 1){
	                int v = e[i].v;
	                cur[u] = i;
	                pre[v] = u;
	                flag = true;
	                u = v;
	                break;
	            }
	        }
	        if(flag){
	            if(u == t){
	                int det = INF;
	                for(int j = u ; j != s ; j = pre[j]){
	                    det = min(det , e[cur[pre[j]]].c);
	                }
	                for(int j = u ; j != s ; j = pre[j]){
	                    e[cur[pre[j]]].c -= det;
	                    e[cur[pre[j]] ^ 1].c += det;
	                }
	                maxflow += det;
	                u = s;
	            }
	 
	        }
	        else{
	            int mind = n;
	            for(int i = head[u] ; i != -1 ; i = e[i].next){
	                if(e[i].c > 0 && dis[e[i].v] < mind){
	                    mind = dis[e[i].v];
	                    cur[u] = i;
	                }
	 
	            }
	            if((--gap[dis[u]]) == 0)break;
	            gap[dis[u] = mind + 1] ++;
	            if(u != s)u = pre[u];
	        }
	    }
	    return maxflow;
	}

#计算几何
##二维

	///考虑误差的加法运算
	double add(double a , double b)
	{
	    if(abs(a + b) < eps * (abs(a) + abs(b)))return 0;
	    return a + b;
	}
	///二维向量结构体
	struct P
	{
	    double x , y;
	    P(){}
	    P(double x , double y) : x(x) , y(y){}
	    P operator + (P p)
	    {
	        return P(add(x , p.x) , add(y , p.y));
	    }
	    P operator - (P p)
	    {
	        return P(add(x , - p.x) , add(y , - p.y));
	    }
	    P operator * (double d)
	    {
	        return P(x * d , y * d);
	    }
	    double dot(P p)///点乘
	    {
	        return add(x * p.x  , y * p.y);
	    }
	    double det(P p)///叉乘
	    {
	        return add(x * p.y , - y * p.x);
	    }
	};
	///判断点q是否在线段p1-p2上
	bool on_seg(P p1 , P p2 , P q)
	{
	    return (p1 - q).det(p2 - q) == 0 && (p1 - q).dot(p2 - q) <= 0;
	}
	///计算直线p1-p2与直线q1-q2的交点坐标
	P intersection(P p1 , P p2 , P q1 , P q2)
	{
	    return p1 + (p2 - p1) * ((q2 - q1).det(q1 - p1) / (q2 - q1).det(p2 - p1));
	}
	///字典序比较
	bool cmp_x(const P& p , const P& q)
	{
	    if(p.x != q.x)return p.x < q.x;
	    return p.y < q.y;
	}
	///求凸包
	vector<P> convex_hull(P* ps , int n)
	{
	    sort(ps , ps + n , cmp_x);
	    int k = 0;///凸包顶点数
	    vector<P> qs(n * 2);///构造中的凸包
	    ///构造凸包的下侧
	    for(int i = 0 ; i < n ; i++)
	    {
	        while(k > 1 && (qs[k - 1] - qs[k - 2]).det(ps[i] - qs[k - 1]) <= 0)k--;
	        qs[k++] = ps[i];
	    }
	    ///构造凸包的上侧
	    for(int i = n - 2 , t = k ; i >= 0 ; i--)
	    {
	        while(k > t && (qs[k - 1] - qs[k - 2]).det(ps[i] - qs[k - 1]) <= 0)k--;
	        qs[k++] = ps[i];
	    }
	    qs.resize(k - 1);
	    return qs;
	}
	///距离的平方
	double dist(P p , P q)
	{
	    return (p - q).dot(p - q);
	}
	int n;
	P ps[MAXN];
	///旋转卡壳
	void solve()
	{
	    vector<P> qs = convex_hull(ps , n);
	    int n = qs.size();
	    if(n == 2)///特别处理凸包退化的情况
	    {
	        printf("%.0f\n" , dist(qs[0] , qs[1]));
	        return ;
	    }
	    int i = 0, j = 0;
	    for(int k = 0  ; k < n ; k++)
	    {
	        if(!cmp_x(qs[i] , qs[k]))i = k;///最右（上）点的序号
	        if(cmp_x(qs[j] , qs[k]))j = k;///最左（下）点的序号
	    }
	    double res = 0;
	    int si = i , sj = j;
	    while(i != sj || j != si)///直到旋转180度
	    {
	        res = max(res , dist(qs[i] , qs[j]));
	        ///判断先转到边i——(i + 1)的法线方向还是边j——(j + 1)的法线方向
	        if((qs[(i + 1) % n] - qs[i]).det(qs[(j + 1) % n] - qs[j]) < 0)
	        {
	            i = (i + 1) % n;
	        }
	        else j = (j + 1) % n;
	    }
	    printf("%.0f\n" , res);
	}
