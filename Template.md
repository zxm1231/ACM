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

###K短路
	//以终点t为源对反图spfa
	struct Edge
	{
	    int v , c , next;
	}e[MAXE] , re[MAXE];
	struct node
	{
	    int now,g,f;
	    bool operator <(const node a)const
	    {
	        if(a.f == f) return a.g < g;
	        return a.f < f;
	    }
	};

	int dis[MAXN];
	int in[MAXN];
	int head[MAXN];
	int rehead[MAXN];
	int cnt[MAXN];
	int num;
	void init()
	{
	    num = 0;
	    clr(head , -1);
	    clr(rehead , -1);
	}
	void adde(int u , int v , int c)
	{
	    e[num].v = v;
	    e[num].c = c;
	    e[num].next = head[u];
	    head[u] = num;
	    re[num].v = u;
	    re[num].c = c;
	    re[num].next = rehead[v];
	    rehead[v] = num++;
	}
	int Astar(int src,int to)
	{
	    priority_queue<node> Q;
	    int i,cnt = 0;
	    if(src == to) k++;//在起点与终点是同一点的情况下，k要+1
	    if(dis[src] == inf) return -1;
	    node a,next;
	    a.now = src;
	    a.g = 0;
	    a.f = dis[src];
	    Q.push(a);
	    while(!Q.empty())
	    {
	        a = Q.top();
	        Q.pop();
	        if(a.now == to)
	        {
	            cnt++;
	            if(cnt == k)
	                return a.g;
	        }
	        for(i = head[a.now]; i!=-1; i = e[i].next)
	        {
	            next = a;
	            next.now = e[i].v;
	            next.g = a.g+e[i].c;
	            next.f = next.g+dis[next.now];
	            Q.push(next);
	
	        }
	    }
	    return -1;//不能到达
	}

##最小生成树
###Prim
	struct Node
	{
	    int v , dis;
	    bool operator < (const Node & n)const{
	        return dis > n.dis;
	    }
	};
	int vis[MAXN];
	int m[MAXN][MAXN]; //邻接矩阵
	int dis[MAXN];
	priority_queue<Node> q;
	int n;
	int prim()
	{
	    int size = 1 , sum = 0;
	    vis[1] = 1;
	    Node node;
	    for(int i = 2 ; i <= n ; i++){
	        dis[i] = m[1][i];
	        node.v = i;
	        node.dis = m[1][i];
	        q.push(node);
	    }
	    while(size < n){
	        node = q.top();
	        q.pop();
	        if(vis[node.v])continue;
	        vis[node.v] = 1;
	        size ++;
	        sum += node.dis;
	        for(int i = 2 ; i <= n ; i++){
	            if(!vis[i] && dis[i] > m[node.v][i]){
	                dis[i] = m[node.v][i];
	                Node tmp;
	                tmp.v = i;
	                tmp.dis = m[node.v][i];
	                q.push(tmp);
	            }
	        }
	    }
	    return sum;
	}

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

##2-SAT
综述：每个条件的形式都是x[i]为真/假或者x[j]为真/假，
每个x[i]拆成2*i和2*i+1两个点，分别表示x[i]为真，x[i]
为假；加的每一条边之间的关系是and
 
模型一：两者（A，B）不能同时取（但可以两个都不选）
说明：A 为假或 B 为假
那么选择了 A 就只能选择 B’，选择了 B 就只能选择 A’
连边 A→B’，B→A’
 
模型二：两者（A，B）不能同时不取（但可以两个都选）
说明：A 为真或 B 为真
那么选择了 A’就只能选择 B，选择了 B’就只能选择 A
连边 A’→B，B’→A
 
模型三：两者（A，B）要么都取，要么都不取
说明：......
那么选择了 A，就只能选择 B，选择了 B 就只能选择 A，选择了 A’就只能选择 B’，
选择了 B’就只能选择 A’
连边 A→B，B→A，A’→B’，B’→A’
 
模型四：两者（A，A’）必取A
那么，那么，该怎么说呢？先说连边吧。
连边 A’→A
 
模型五（补充） ：两者（A，B）两个必须不相同，即要么选A，要么选B
逻辑表达：A||B 非 A||非 B
连边：A 为真或 B 为真： A’--->B B’--->A;
A 为假或 B 为假: A-->B’ B-->A
说明：A 或 B，非 A 或非 B，前者表示两者至少有一个 true，后者表示至少有一个 false

###2-SAT + 二分答案  统一建模的方式：

同一组的两个状态分别存储在2*i和2*i+1两个节点，产生2*n个节点

	for(int i=1;i<2*n;i++)
		for(int j=0;j<i;j++)
		{
			if (i==(j^1)) continue;//记得j^1加上小括号
			sat.add_clause(i,j);//枚举出的不属于同一组的不相容的两点
		}
###【O(NM)算法：求字典序最小的解】
根据2-SAT建成的图中边的定义可以发现，若图中i到j有路径，则若i选，则j也要选；或者说，若j不选，则i也不能选；
因此得到一个很直观的算法：
（1）给每个点设置一个状态V，V=0表示未确定，V=1表示确定选取，V=2表示确定不选取。称一个点是已确定的当且仅当其V值非0。设立两个队列Q1和Q2，分别存放本次尝试选取的点的编号和尝试不选的点的编号。
（2）若图中所有的点均已确定，则找到一组解，结束，否则，将Q1、Q2清空，并任选一个未确定的点i，将i加入队列Q1，将i'加入队列Q2；
（3）找到i的所有后继。对于后继j，若j未确定，则将j加入队列Q1；若j'（这里的j'是指与j在同一对的另一个点）未确定，则将j'加入队列Q2；
（4）遍历Q2中的每个点，找到该点的所有前趋（这里需要先建一个补图），若该前趋未确定，则将其加入队列Q2；
（5）在（3）（4）步操作中，出现以下情况之一，则本次尝试失败，否则本次尝试成功：
<1>某个已被加入队列Q1的点被加入队列Q2；
<2>某个已被加入队列Q2的点被加入队列Q1;
<3>某个j的状态为2；
<4>某个i'或j'的状态为1或某个i'或j'的前趋的状态为1；
（6）若本次尝试成功，则将Q1中的所有点的状态改为1，将Q2中所有点的状态改为2，转（2），否则尝试点i'，若仍失败则问题无解。
该算法的时间复杂度为O(NM)（最坏情况下要尝试所有的点，每次尝试要遍历所有的边），但是在多数情况下，远远达不到这个上界。
具体实现时，可以用一个数组vst来表示队列Q1和Q2。设立两个标志变量i1和i2（要求对于不同的i，i1和i2均不同，这样可以避免每次尝试都要初始化一次，节省时间），若vst[i]=i1则表示i已被加入Q1，若vst[i]=i2则表示i已被加入Q2。不过Q1和Q2仍然是要设立的，因为遍历（BFS）的时候需要队列，为了防止重复遍历，加入Q1（或Q2）中的点的vst值必然不等于i1（或i2）。中间一旦发生矛盾，立即中止尝试，宣告失败。

该算法虽然在多数情况下时间复杂度到不了O(NM)，但是综合性能仍然不如下面的O(M)算法。不过，该算法有一个很重要的用处：求字典序最小的解！
如果原图中的同一对点编号都是连续的（01、23、45……）则可以依次尝试第0对、第1对……点，每对点中先尝试编号小的，若失败再尝试编号大的。这样一定能求出字典序最小的解（如果有解的话），因为一个点一旦被确定，则不可更改。
如果原图中的同一对点编号不连续（比如03、25、14……）则按照该对点中编号小的点的编号递增顺序将每对点排序，然后依次扫描排序后的每对点，先尝试其编号小的点，若成功则将这个点选上，否则尝试编号大的点，若成功则选上，否则（都失败）无解。

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
##最小费用流
###spfa增广费用流
	struct Edge
	{
	    int from;
	    int to;
	    int next;
	    int re;//记录逆边的下标
	    int cap;//容量
	    int cost;//费用
	}e[MAXE];
	int pre[MAXN];
	int head[MAXN];
	bool in[MAXN];
	int que[MAXN];
	int dis[MAXN];
	int num;//边的总数
	void init()
	{
	    num = 0;
	    clr(head , -1);
	}
	void add(int u,int v,int ca,int co)
	{
	    e[num].from=u;
	    e[num].to=v;
	    e[num].cap=ca;
	    e[num].cost=co;
	    e[num].re=num+1;
	    e[num].next=head[u];
	    head[u]=num++;
	 
	    e[num].from=v;//加逆边
	    e[num].to=u;
	    e[num].cap=0;
	    e[num].cost=-co;
	    e[num].re=num-1;
	    e[num].next=head[v];
	    head[v]=num++;
	}
	int n;
	int start;
	int end;
	bool SPFA()
	{
	    int front=0,rear=0;
	    for(int v=0;v<=n;v++)
	    {
	        if(v==start)
	        {
	            que[rear++]=v;
	            in[v]=true;
	            dis[v]=0;
	        }
	        else
	        {
	            dis[v]=INF;
	            in[v]=false;
	        }
	    }
	    while(front!=rear)
	    {
	        int u=que[front++];
	        in[u]=false;
	        if(front>=MAXN)front=0;
	        for(int i=head[u];i!=-1;i=e[i].next)
	        {
	            int v=e[i].to;
	            if(e[i].cap&&dis[v]>dis[u]+e[i].cost)
	            {
	                dis[v]=dis[u]+e[i].cost;
	                pre[v]=i;
	                if(!in[v])
	                {
	                    que[rear++]=v;
	                    in[v]=true;
	                    if(rear>=MAXN)rear=0;
	                }
	            }
	        }
	    }
	    if(dis[end]==INF)return false;
	    return true;
	}
	int c;//费用
	int f;//最大流
	 
	void minCostMaxflow()
	{
	    c=f=0;
	    int u,p;
	    while(SPFA())
	    {
	        int Min=INF;
	        for(u=end;u!=start;u=e[p].from)
	        {
	            p=pre[u];
	            Min=min(Min,e[p].cap);
	        }
	        for(u=end;u!=start;u=e[p].from)
	        {
	            p=pre[u];
	            e[p].cap-=Min;
	            e[e[p].re].cap+=Min;
	 
	        }
	        c+=dis[end]*Min;
	        f+=Min;
	    }
	}