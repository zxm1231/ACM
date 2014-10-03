#配置环境

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


#图论
##最短路
###Dijstra
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
	int lowcost[MAXN],closest[MAXN];
	int prim(int v0)
	{
	    int i,j,mindis,minone;
	    int ans = 0;/*用来记录最小生成树的总长度*/
	    /*各点距离初始化*/
	    for(i = 0;i < n;i++)
	    {
	        lowcost[i] = cost[v0][i];
	        closest[i] = v0;
	    }
	    for(i = 0;i < n-1;i++)
	    {
	        mindis = inf;
	        for(j = 0;j < n;j++)
	          if(lowcost[j] && mindis > lowcost[j])
	          {
	              mindis = lowcost[j];
	              minone = j;
	          }
	        /*将找到的最近点加入最小生成树*/
	        ans += lowcost[minone];
	        lowcost[minone] = 0;
	        /*修正其他点到最小生成树的距离*/
	        for(j = 0;j < n;j++)
	          if(cost[j][minone] < lowcost[j])
	          {
	              lowcost[j] = cost[j][minone];//邻接矩阵
	              closest[j] = minone;
	          }
	    }
	    return ans;
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
##拓扑排序
###1.朴素算法：
将入度为0的点加入队列，再依次取出队列中的元素，把其出边依次遍历，并且将边指向的节点入度减一，同时将入度为0的点加入队列。

	queue < int > Q;
	memset(idegree, 0, sizeof(idegree));
	memset(last, -1, sizeof(last));
	tot = 0;
	for(int i = 0; i < m; i++)
	{
		int x, y;
		scanf("%d%d", &x, &y);
		addedge(x, y);
		idegree[y]++;
	}
	int now = 0;
	for(int i = 1; i <= n; i++)
		if(idegree[i] == 0)
			Q.push(i);
	while(!Q.empty())
	{
		int u = Q.top();
		Q.pop();
		now++;
		printf("%d", u);
		if(now < n) putchar(' ');
		else putchar('\n');
		for(int j = last[u]; -1 != j; j = edge[j].next)
		{
			int v = edge[j].v;
			idegree[v]--;
			if(idegree[v] == 0)
				Q.push(v);
		}
	}

###2.求字典序最小的方案：

		当各点标号均互不相同时，直接将队列改为优先队列；
		当存在重复标号时？？？？
		
###3.基于DFS的拓扑排序（很短很好写）
		记录各点完成访问的时刻(完成时间),用DFS遍历一次整个图,得出各结点的完成时间,然后按完成时间倒序排列就得到了图的拓扑序列

	/* 拓扑排序         O(e)
	 * 确保是有向无环图！
	 * 结果逆序存放在sta中！
	 * */
	VI ve[Maxn];
	bool vis[Maxn];
	int sta[Maxn], top;
	void dfsTopo(int u) {
	    vis[u] = true;
	    for (int i = 0; i < ve[u].size(); i ++ )
	           if (!vis[ve[u][i]])
	               dfsTopo(ve[u][i]);
	    sta[top ++ ] = u;
	}
	void Toposort(int n) {
	    memset(vis, 0, sizeof(vis));
	    top = 0;
	    for (int i = 1; i <= n; i ++ )
	        if (!vis[i]) dfsTopo(i);
	}
###4.拓扑排序数
求树的拓扑排序数：
dp[root] = num[root] ! / (num[i] , i in tree[root])

###5.关键路径
求关键路径：
#### 关键节点：
		正向拓扑序，求每个点最早到达时间early[i] = max(early[i], early[j] + edge[k]);
		利用汇点的early值，反向拓扑求每个点迟到达时间late[i] = min(late[i], late[j] - edge[k]);
		关键节点early == late
		PS：
				最早开始时间可以理解为是必须等到之前的任务完成才能做
				最迟开始时间可以理解为是必须为后面的任务留出足够的时间
		
#### 关键路径：
	源点到汇点的最长路（可能有多条）
	/*
	* 求出early 和 latest
	* idegree  odegree 出入度
	* 中间有重建图，所以保存了一条边的u 和 v
	*/
	
	memset(idegree, 0, sizeof(idegree));
	memset(odegree, 0, sizeof(odegree));
	memset(last, -1, sizeof(last));
	tot = 0;
	for(int i = 0; i < m; i++)
	{
		int x, y, z;
		scanf("%d%d%d", &x, &y, &z);
		addedge(x, y, z);
		idegree[y]++;
		odegree[x]++;
	}
	memset(early, 0, sizeof(early));
	for(int i = 1; i <= n; i++)
		if(!idegree[i])
			Q.push(i);
	while(!Q.empty())
	{
		int u = Q.front();
		Q.pop();
		for(int j = last[u]; -1 != j; j = edge[j].next)
		{
			int v = edge[j].v;
			idegree[v]--;
			if(!idegree[v])
				Q.push(v);
			int temp = early[u] + edge[j].len;
			early[v] = max(early[v], temp);
		}
	}
	memset(last, -1, sizeof(last));
	tot = 0;
	for(int i = 0; i < m; i++)
		addedge(edge[i].v, edge[i].u, edge[i].len);
	memset(latest, CLRINF, sizeof(latest));
	cout << latest[0] << endl;
	for(int i = 1; i <= n; i++)
		if(odegree[i] == 0)
		{
			latest[i] = early[i];
			Q.push(i);
		}
	while(!Q.empty())
	{
		int u = Q.front();
		Q.pop();
		for(int j = last[u]; -1 != j; j = edge[j].next)
		{
			int v = edge[j].v;
			odegree[v]--;
			if(!odegree[v])
				Q.push(v);
			int temp = latest[u] - edge[j].len;
			latest[v] = min(latest[v], temp);
		}
	}

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

###点双联通分量
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

###无向图全局最小割

	int n , m;
	int combine[MAXN];
	int g[MAXN][MAXN] , node[MAXN];
	int st , en , minCut, k;
	int top, sta[MAXN];
	int maxi;
	int vis[MAXN];
	int wet[MAXN];
	int Search(int n)
	{
	    clr(vis,0);
	    clr(wet,0);
	    int minCut = 0;
	    int temp = -1;
	    st = -1, en = -1;
	    int top = 0;
	    for(int i=0; i< n; i++)
	    {
	        int maxi = -INF;
	        for(int j = 0; j < n; j++)
	        {
	            int u = node[j];
	            if(!combine[u] && !vis[u] && wet[u] > maxi)
	            {
	                temp = u;
	                maxi = wet[u];
	            }
	        }
	        sta[top++] = temp;
	        vis[temp] = true;
	        if(i == n - 1)
	            minCut = maxi;
	        for(int j = 0; j < n; j++)
	        {
	            int u = node[j];
	            if(!combine[u] && !vis[u])
	            {
	                wet[u] += g[temp][u];
	            }
	        }
	    }
	    st = sta[top - 2];
	    en = sta[top - 1];
	    for(int i = 0; i < top; i++)  node[i] = sta[i];
	    return minCut;
	}
	
	int SW(int n)
	{
	    int ans = inf;
	    clr(combine,0);
	    for(int i = 0; i < n; i++)
	        node[i] = i;
	    for(int i = 1; i < n; i++)
	    {
	        k = n - i + 1;
	        int cur = Search(k);
	        if(cur < ans)
	        {
	            ans = cur;
	        }
	        if(ans == 0) return ans;
	        combine[en] = true;
	        for(int j = 0; j < n; j++)
	        {
	            if(j == st) continue;
	            if(!combine[j])
	            {
	                g[st][j] += g[en][j];
	                g[j][st] += g[j][en];
	            }
	        }
	    }
	    return ans;
	}
	
	int main()
	{
	    while(~scanf("%d%d" , &n , &m)){
	        clr(g , 0);
	        for(int i = 0 ; i < m ; i++){
	            int u , v , c;
	            scanf("%d%d%d" , &u , &v , &c);
	            g[u][v] += c;
	            g[v][u] += c;
	        }
	        printf("%d\n" , SW(n));
	    }
	    return 0;
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
###一般图最大匹配——带花树
	#define N 1100
	#define M 200100
	#define LEN 100100
	#define INF (1 << 30)
	typedef long long LL;
	
	int n, head, tail, Start, Finish;
	int match[N];     //表示哪个点匹配了哪个点
	int Father[N];   //这个就是增广路的father……但是用起来太精髓了
	int base[N];     //该点属于哪朵花
	int Q[N];
	bool mark[N], map[N][N], InBlossom[N], in_Queue[N];
	
	void init()
	{
	    int x,y;
	    scanf("%d",&n);
	    while (scanf("%d%d",&x,&y)!=EOF)
	        map[x][y]=map[y][x]=1;
	}
	
	void BlossomContract(int x,int y)
	{
	    clr(mark, false);
	    clr(InBlossom,false);
	#define pre father[match[i]]
	    int lca,i;
	    for (i=x; i; i=pre)
	    {
	        i=base[i];
	        mark[i]=true;
	    }
	    for (i=y; i; i=pre)
	    {
	        i=base[i];    //寻找lca之旅……一定要注意i=base[i]
	        if (mark[i])
	        {
	            lca=i;
	            break;
	        }
	    }
	    for (i=x; base[i]!=lca; i=pre)
	    {
	        if (base[pre]!=lca) father[pre]=match[i]; //对于BFS树中的父边是匹配边的点，father向后跳
	        InBlossom[base[i]]=true;
	        InBlossom[base[match[i]]]=true;
	    }
	    for (i=y; base[i]!=lca; i=pre)
	    {
	        if (base[pre]!=lca) father[pre]=match[i]; //同理
	        InBlossom[base[i]]=true;
	        InBlossom[base[match[i]]]=true;
	    }
	#undef pre
	    if (base[x]!=lca) father[x]=y;     //注意不能从lca这个奇环的关键点跳回来
	    if (base[y]!=lca) father[y]=x;
	    for (i=1; i<=n; i++)
	        if (InBlossom[base[i]])
	        {
	            base[i]=lca;
	            if (!in_Queue[i])
	            {
	                Q[++tail]=i;
	                in_Queue[i]=true;     //要注意如果本来连向BFS树中父结点的边是非匹配边的点，可能是没有入队的
	            }
	        }
	}
	
	void Change()
	{
	    int x,y,z;
	    z=Finish;
	    while (z)
	    {
	        y=father[z];
	        x=match[y];
	        match[y]=z;
	        match[z]=y;
	        z=x;
	    }
	}
	
	void FindAugmentPath()
	{
	    clr(father, 0);
	    clr(in_Queue, false);
	    for (int i=1; i<=n; i++) base[i]=i;
	    head=0;
	    tail=1;
	    Q[1]=Start;
	    in_Queue[Start]=1;
	    while (head!=tail)
	    {
	        int x=Q[++head];
	        for (int y=1; y<=n; y++)
	            if (map[x][y] && base[x]!=base[y] && match[x]!=y)   //无意义的边
	                if ( Start==y || match[y] && father[match[y]] )    //精髓地用father表示该点是否
	                    BlossomContract(x,y);
	                else if (!father[y])
	                {
	                    father[y]=x;
	                    if (match[y])
	                    {
	                        Q[++tail]=match[y];
	                        in_Queue[match[y]]=true;
	                    }
	                    else
	                    {
	                        Finish=y;
	                        Change();
	                        return;
	                    }
	                }
	    }
	}
	
	void Edmonds()
	{
	    clr(match, 0);
	    for (Start=1; Start<=n; Start++)
	        if (match[Start]==0)
	            FindAugmentPath();
	}
	
	void output()
	{
	    clr(mark, false);
	    int cnt=0;
	    for (int i=1; i<=n; i++)
	        if (match[i]) cnt++;
	    printf("%d\n",cnt);
	    for (int i=1; i<=n; i++)
	        if (!mark[i] && match[i])
	        {
	            mark[i]=true;
	            mark[match[i]]=true;
	            printf("%d %d\n",i,match[i]);
	        }
	}
	
	int main()
	{
	//    freopen("input.txt","r",stdin);
	    init();
	    Edmonds();
	    output();
	    return 0;
	}

##欧拉回路
	/*欧拉回路, 有向图*/
	vec ve[Maxn];
	int cur[Maxn];
	stack<int> eulerianWalk(int u) { //返回欧拉回路的逆序
	    stack<int> sta, ret;
	    sta.push(u);
	    cur[u] = 0;
	    while (!sta.empty()) {
	        u = sta.top();
	        sta.pop();
	        while (cur[u] < ve[u].size()) {
	            sta.push(u);
	            u = ve[u][cur[u] ++ ];
	        }
	        ret.push(u);
	    }
	    return ret;
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
###ZKW流
	struct node
	{
	    int u, v, c, w, next;
	}edge[M];
	int tot, last[N];
	int dist[N],pre[N];
	bool visit[N];
	int n, m, src, des;
	int flow, cost, value;
	
	void addedge(int u, int v, int c, int w)
	{
	    edge[tot].u = u; edge[tot].v = v; edge[tot].c = c; edge[tot].w = w; edge[tot].next = last[u]; last[u] = tot++;
	    edge[tot].u = v; edge[tot].v = u; edge[tot].c = 0; edge[tot].w = -w; edge[tot].next = last[v]; last[v] = tot++;
	}
	
	int Aug(int u, int m)
	{
		if(u == des)
		{
			cost += value * m;
			flow += m;
			return m;
		}
		visit[u] = true;
		int l = m;
		for(int j = last[u]; j != -1; j = edge[j].next)
		{
		    int v = edge[j].v, c = edge[j].c, w = edge[j].w;
			if(c && !w && !visit[v])
			{
				int delta = Aug(v, l < c ? l : c);
				edge[j].c -= delta;
				edge[j ^ 1].c += delta;
				l -= delta;
				if(!l) return m;
			}
		}
		return(m - l);
	}
	
	bool ModLabel(int src, int des)
	{
		int d = INF;
		for(int i = 0; i <= des; i++)
			if(visit[i])
			{
				for(int j = last[i]; j != -1; j = edge[j].next)
				{
					if(edge[j].c && !visit[edge[j].v] && edge[j].w < d) d = edge[j].w;
				}
			}
		if(d == INF) return false;
		for(int i = 0; i <= des; i++)
		  if(visit[i])
		  {
			  for(int j = last[i]; j != -1; j = edge[j].next)
			  {
				  edge[j].w -= d;
				  edge[j^1].w += d;
			   }
		  }
		value += d;
		return true;
	}
	
	void MinCostMaxFlow(int src, int des)
	{
	    flow = cost = value = 0;
	    int xx = 0;
	    do
	    {
	        //cout << xx ++ << endl;
	        do
	        {
	            memset(visit, 0, sizeof(visit));
	        }while(Aug(src, INF));
	    }while(ModLabel(src, des));
	}
###注意：
B[u,v]表示(u,v)流量的下限，C[u,v]表示(u,v)流量的上限, F[u,v]表示(u,v)的流量, 
g[u,v]表示F[u,v]-B[u,v] 显然 0<=g[u,v]<=C[u,v]-B[u,v]  

####1.无源汇的可行流 :
我们要想办法转换为有源汇的最大流问题.
考虑流量都为g[,]且容量为C[,]-B[,]的网络，貌似有点接近最后的转换方式了，
为了不忽略B[,]这一条件，我们把g[,]最后强制加上B[,].
但会发现一个致命漏洞，加上后就未必满足流量平衡了！
对于这个有两种办法解决..
一种方法是添加附加源汇S,T  对于某点 u, 设 M(u)=sigma(B[i,u])-sigma(B[u,j]) ，
则根据流量平衡条件有 M(u)同时等于 sigma(g[u,j])-sigma(g[i,u])
若M(u)<0，即sigma(g[u,j]) < sigma(g[i,u]) 进入u的流量比从u 出去的多，
所以 u -> T 连容量为  -(sigma(B[i,u])-sigma(B[u,j]) ) 的边
同理. M(u)>0时，即 S->u 连容量为 sigma(B[i,u])-sigma(B[u,j])  的边.
然后再对于任意边(i,u)/(u,j) 连一条 C[u,v]-B[u,v]的边.
这样只需对新的网络求一遍最大流即可. 若出附加源点的边都满流即是存在可行流，反之不然.
满流的必要条件是显然的. 不满流不能保证加上B[,]后流量平衡. 前面都白费了.
另一种方法相对简单.其实类似，本质相同.
仍添加附加源汇S,T 对于某边 (u,v) 在新网络中连边 S->v 容量 B[u,v]   ,  u->T 容量 B[u,v]  ， u->v 容量 C[u,v]-B[u,v]
可以这样理解，边S->v : 求的时候直接从S流过来的流量值B[u,v], 与最终解中边(u,v)强制加上的从 u流过来的流量B[u,v]，对v点的流量平衡条件的影响实质等价.
边u->T同理.
最后，一样也是求一下新网络的最大流，判断从附加源点的边，是否都满流即可.

具体的解？根据最前面提出的强制转换方式，边(u,v)的最终解中的实际流量即为g[u,v]+B[u,v]
为什么这种方法只适用于无源汇上下界可行流？
本质上是因为S,T并不满足流量平衡，而上述的方法都是考虑到每点的流量平衡而建的. 但有些时候貌似还是可以出正确解. 至于有没有什么解决方法，下次再想想吧~【标记下】
例题 ZOJ 2314 / SGU 194 Reactor Cooling http://acm.sgu.ru/problem.php?contest=0&problem=194 
####2.有源汇的上下界可行流
从汇点到源点连一条上限为INF，下限为0的边. 按照 1.无源汇的上下界可行流一样做即可.
改成无源汇后，求的可行流是类似环的，流量即T->S边上的流量.  这样做显然使S,T也变得流量平衡了.
####3.有源汇的上下界最大流
方法一： 2.有源汇上下界可行流中，从汇点到源点的边改为连一条上限为INF，下限为x的边.
因为显然x>ans即MIN(T->S )> MAX(S->T) ,会使求新网络的无源汇可行流无解的（S,T流量怎样都不能平衡）
而x<=ans会有解.
所以满足二分性质，二分x，最大的x使得新网络有解的即是所求答案原图最大流.
方法二：从汇点T到源点S连一条上限为INF，下限为0的边，变成无源汇的网络.  照求无源汇可行流的方法(如1)，建附加源点S'与汇点T'，求一遍S'->T‘的最大流. 再把从汇点T到源点S的这条边拆掉 . 求一次从S 到T 的最大流即可. （关于S',T'的边好像可以不拆？）（这样一定满足流量平衡？）表示这方法我也没有怎么理解.
####4.有源汇的上下界最小流
方法一： 2.有源汇上下界可行流中，从汇点到源点的边改为连一条上限为x，下限为0的边.
与3同理，二分上限，最小的x使新网络无源汇可行流有解，即是所求答案原图最小流.
方法二:  照求无源汇可行流的方法(如1)，建附加源点S'与汇点T'，求一遍S'->T‘的最大流. 但是注意这一遍不加汇点到源点的这条边，即不使之改为无源汇的网络去求解. 求完后，再加上那条汇点到源点上限INF的边. 因为这条边下限为0，所以S',T'无影响. 再直接求一遍S'->T'的最大流. 若S’出去的边全满流，T->S边上的流量即为答案原图最小流，否则若不全满流即无解. 
和求3.有源汇的上下界最大流过程相反，感性理解是:  
首先明确，我们的方法是通过加边转化成对任一点都有流量平衡的无源汇的网络，进行求解.
即最终解只能是加上边后，求的无源汇可行流，即T->S这边上的流量.  不改成无源汇的直接求的解是未必正确的，在（1）中已经提到.
然后，因为第一遍做的时候并无这条边，所以S->T的流量在第一遍做的时候都已经尽力往其他边流了. 于是加上T->S这条边后，都是些剩余的流不到其他边的流量. 从而达到尽可能减少T->S这边上的流量的效果，即减小了最终答案.
感觉上第一遍做的既然是不改成无源汇直接求的，应该是错误的？
这里不是错误的. 首先我们的解都是按照第二遍所求的而定，其次这里这样做本质是延迟对T->S这条边的增流.
