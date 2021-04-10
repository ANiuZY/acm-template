# ANiu 的ACM模板


## 数据结构



```c++
//__builtin_popcount(); 二进制数位1的个数。
//值域分块：支持O(1)查询前缀和，不超过O(sqrt(n))的时间进行单点修改。
```



### 快速区间操作 

```c++
//log修改 log查询。
struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int val){
        for(;x<=n;x+=lowbit(x))a[x]+=val;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))res+=a[x];
        return res;
    }
}tr;

//sqrt(n)修改，o(1)查询。
struct BT{
    int k[M],as[M],aq[M],s,ns;
    void init(){// 使用前先init()!!!
        s=sqrt(n);ns=n/s+1;
        for(int i=1;i<=n;++i)k[i]=i/s+1;
    }
    void add(int x,int val){
        for(int i=k[x]+1;i<=ns;++i)as[i]+=val;
        for(int x0=x;x<=n&&k[x]==k[x0];++x)aq[x]+=val;
    }
    int get(int x){
        return as[k[x]]+aq[x];
    }
}bt;
```



### RMQ

```c++
#include<bits/stdc++.h>
#define M 10000
using namespace std;
int d[M][10]

void RMQ_init(const vector<int>& A){
    int n=A.size();
    for(int i=0;i<n;++i)d[i][0]=A[i];
    for(int j=1;(1<<j)<=n;++j)
        for(int i=0;i+(1<<j)<=n;++i)
            d[i][j]=min(d[i][j-1],d[i+(1<<(j-1))][j-1]);
}
int RMQ(int L,int R){
    int k=0;
    while((1<<(k+1))<=R-L+1)++k;
    return min(d[L][k],d[R-(1<<k)+1][k]);
}
```



### 线段树 - 区间和、积

```c++
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
typedef long long ll;
ll n,m,p;
ll a[M];

struct nd{
    ll sum,lzs,lzm;
    ll len;
}tr[M<<2];
ll L,R,val;

void pu(ll o){
    tr[o].sum=(tr[o<<1].sum+tr[o<<1|1].sum)%p;
}

void pd(ll o){
    if(tr[o].lzm!=1){
        ll v=tr[o].lzm;
        tr[o<<1].sum=(tr[o<<1].sum*v)%p;
        tr[o<<1|1].sum=(tr[o<<1|1].sum*v)%p;
        tr[o<<1].lzm=(tr[o<<1].lzm*v)%p;
        tr[o<<1].lzs=(tr[o<<1].lzs*v)%p;
        tr[o<<1|1].lzm=(tr[o<<1|1].lzm*v)%p;
        tr[o<<1|1].lzs=(tr[o<<1|1].lzs*v)%p;
        tr[o].lzm=1;
    }
    if(tr[o].lzs){
        ll v=tr[o].lzs;
        tr[o<<1].sum=(tr[o<<1].sum+v*tr[o<<1].len)%p;
        tr[o<<1|1].sum=(tr[o<<1|1].sum+v*tr[o<<1|1].len)%p;
        tr[o<<1].lzs=(tr[o<<1].lzs+v)%p;
        tr[o<<1|1].lzs=(tr[o<<1|1].lzs+v)%p;
        tr[o].lzs=0;
    }
}

void build(ll o,ll l,ll r){
    tr[o].lzm=1;
    if(l==r){
        tr[o].sum=a[l];
        tr[o].len=1;
        return;
    }
    ll m=l+r>>1;
    build(o<<1,l,m);
    build(o<<1|1,m+1,r);
    tr[o].len=tr[o<<1].len+tr[o<<1|1].len;
    pu(o);
}

void add(ll o,ll l,ll r){
    if(r<L||l>R)return;
    if(r<=R&&l>=L){
        tr[o].sum=(tr[o].sum+val*(r-l+1))%p;
        tr[o].lzs=(tr[o].lzs+val)%p;
    }else{
        pd(o);
        ll m=l+r>>1;
        add(o<<1,l,m);
        add(o<<1|1,m+1,r);
        pu(o);
    }
}

void mul(ll o,ll l,ll r){
    if(r<L||l>R)return;
    if(r<=R&&l>=L){
        tr[o].sum=tr[o].sum*val%p;
        tr[o].lzm=tr[o].lzm*val%p;
        tr[o].lzs=tr[o].lzs*val%p;
    }else{
        pd(o);
        ll m=l+r>>1;
        mul(o<<1,l,m);
        mul(o<<1|1,m+1,r);
        pu(o);
    }
}

ll get(ll o,ll l,ll r){
    if(l>R||r<L)return 0;
    if(l>=L&&r<=R)return tr[o].sum;
    pd(o);
    ll m=l+r>>1;
    return (get(o<<1,l,m)+get(o<<1|1,m+1,r))%p;
}

int main(){
    scanf("%lld %lld %lld",&n,&m,&p);
    for(int i=1;i<=n;++i)scanf("%lld",a+i);
    build(1,1,n);
    ll x,y,z;
    for(int i=0;i<m;++i){
        scanf("%lld",&x);
        if(x==1){
            scanf("%lld %lld %lld",&x,&y,&z);
            L=x;R=y;val=z;
            mul(1,1,n);
        }else if(x==2){
            scanf("%lld %lld %lld",&x,&y,&z);
            L=x;R=y;val=z;
            add(1,1,n);
        }else{
            scanf("%lld %lld",&x,&y);
            L=x;R=y;
            printf("%lld\n",get(1,1,n));
        }
    }
}
```



### 线段树 - 扫描线

```c++
#include<bits/stdc++.h>
#define ls o<<1
#define rs o<<1|1
using namespace std;
typedef long long ll;
const int M=400010;
int n,rk[M],top,cnt;
struct nd{
    int x,y1,y2,f;
    nd(){}
    nd(int x1,int x2,int x3,int x4){x=x1;y1=x2;y2=x3;f=x4;}
    bool operator < (const nd& t)const {return x<t.x;}
}a[M];
//扫描线线段树
int x[M<<2],num[M<<2],v[M<<2];
void pu(int o){
    if(num[o])x[o]=v[o];
    else x[o]=x[ls]+x[rs];
}
void build(int o,int l,int r){
    if(l==r){
        v[o]=rk[l+1]-rk[l];
    }else{
        int m=l+r>>1;
        build(ls,l,m);
        build(rs,m+1,r);
        v[o]=v[ls]+v[rs];
    }
}
int L,R,val;
void add(int o,int l,int r){
    if(l>R||r<L)return;
    if(l>=L&&r<=R){
        num[o]+=val;
    }else{
        int m=l+r>>1;
        add(ls,l,m);
        add(rs,m+1,r);
    }
    pu(o);// !
}

int main(){
    scanf("%d",&n);
    int x1,x2,y1,y2;
    for(int i=1;i<=n;++i){
        scanf("%d %d %d %d",&x1,&y1,&x2,&y2);
        rk[++top]=y1;
        rk[++top]=y2;
        a[++cnt]=nd(x1,y1,y2,1);
        a[++cnt]=nd(x2,y1,y2,-1);
    }
    sort(a+1,a+1+cnt);
    sort(rk+1,rk+1+top);
    top=unique(rk+1,rk+1+top)-rk-1;
    build(1,1,top-1);
    ll ans=0;
    int p=1;
    for(int i=1;i<=cnt;++i){
        a[i].y1=lower_bound(rk+1,rk+1+top,a[i].y1)-rk;
        a[i].y2=lower_bound(rk+1,rk+1+top,a[i].y2)-rk;
    }
    for(int i=1;i<=cnt;++i){
        ans+=1ll*(a[i].x-a[i-1].x)*x[1];
        L=a[i].y1;R=a[i].y2-1;val=a[i].f;
        add(1,1,top-1);
    }
    printf("%lld\n",ans);
}
```



### 线段树 - 区间线段覆盖

```c++
//不考虑交叉点;李超线段树;
#include<bits/stdc++.h>
using namespace std;
const int mod=39989; //x轴范围
const int mod1=1e9;  //y轴范围
int n,ans,cnt;
double h;

struct nd{
    int id;
    double y1,y2;
    nd(){}
    nd(double x1,double x2,int x3){y1=x1;y2=x2;id=x3;}
    double get(int x,int l,int r){
        if(l==r)return max(y1,y2);
        return y1+(y2-y1)*(x-l)/(r-l);
    }
}a[mod<<2];

#define ls o<<1
#define rs o<<1|1
nd k;
int L,R,p;
void add(int o,int l,int r){ // 添加一条[L,R]内的线段k
    if(l>R||r<L)return;
    if(l>=L&&r<=R){
        double y1=k.get(l,L,R),y2=k.get(r,L,R);
        if(y1<a[o].y1&&y2<a[o].y2)return;
        if(y1>=a[o].y1&&y2>=a[o].y2){
            a[o]=nd(y1,y2,k.id);
            return;
        }
    }
    int m=l+r>>1;
    add(ls,l,m);
    add(rs,m+1,r);
}
void get(int o,int l,int r){ // x=p处的最大值
    double t=a[o].get(p,l,r);
    if(t>h)ans=a[o].id,h=t;
    if(l==r)return;
    int m=l+r>>1;
    if(p<=m)get(ls,l,m);
    else get(rs,m+1,r);
}

int main(){
    scanf("%d",&n);
    int op,y1,y2;
    while(n--){
        scanf("%d",&op);
        if(op){
            scanf("%d %d %d %d",&L,&y1,&R,&y2);//强制在线
            k.id=++cnt;
            L=(L+ans-1)%mod+1;
            R=(R+ans-1)%mod+1;
            k.y1=(y1+ans-1)%mod1+1;
            k.y2=(y2+ans-1)%mod1+1;
            if(L>R)swap(L,R),swap(k.y1,k.y2);
            add(1,1,mod);
        }else{
            scanf("%d",&p);
            p=(p+ans-1)%mod+1;
            ans=0;h=0;
            get(1,1,mod);
            printf("%d\n",ans);
        }
    }
}
```





### 线段树 - 01序列操作 

```c++
//区间置0置1 区间翻转 区间求和 区间最大连续段
#include<bits/stdc++.h>
#define ls o<<1
#define rs o<<1|1
using namespace std;
const int M=100010;
int n,m,a[M];

int L,R;
struct Seg{
    int one[M<<2],len[M<<2];
    int to0[M<<2],to1[M<<2],rv[M<<2],lone[M<<2],rone[M<<2],mone[M<<2];
    int lzr[M<<2],rzr[M<<2],mzr[M<<2];

    void pu(int o){
        one[o]=one[ls]+one[rs];
        if(one[ls]==len[ls])lone[o]=lone[ls]+lone[rs];
        else lone[o]=lone[ls];
        if(one[rs]==len[rs])rone[o]=rone[rs]+rone[ls];
        else rone[o]=rone[rs];
        mone[o]=max(mone[ls],mone[rs]);
        mone[o]=max(mone[o],rone[ls]+lone[rs]);

        if(one[ls]==0)lzr[o]=lzr[ls]+lzr[rs];
        else lzr[o]=lzr[ls];
        if(one[rs]==0)rzr[o]=rzr[rs]+rzr[ls];
        else rzr[o]=rzr[rs];
        mzr[o]=max(mzr[ls],mzr[rs]);
        mzr[o]=max(mzr[o],rzr[ls]+lzr[rs]);
    }
    void To1(int o){
        to1[o]=1;
        to0[o]=0;
        rv[o]=0;
        lone[o]=rone[o]=mone[o]=one[o]=len[o];
        lzr[o]=rzr[o]=mzr[o]=0;
    }
    void To0(int o){
        to0[o]=1;
        to1[o]=0;
        rv[o]=0;
        lone[o]=rone[o]=mone[o]=one[o]=0;
        lzr[o]=rzr[o]=mzr[o]=len[o];
    }
    void Rs(int o){
        rv[o]++;
        one[o]=len[o]-one[o];
        swap(lone[o],lzr[o]);
        swap(rone[o],rzr[o]);
        swap(mone[o],mzr[o]);
        swap(to1[o],to0[o]);
    }
    void pd(int o){
        if(rv[o]&1){
            rv[o]=0;
            Rs(ls);
            Rs(rs);
        }
        if(to1[o]){
            to1[o]=0;
            To1(ls);
            To1(rs);
        }
        if(to0[o]){
            to0[o]=0;
            To0(ls);
            To0(rs);
        }
    }
    void build(int o,int l,int r){
        if(l==r){
            if(a[l]==1)one[o]=1;
            len[o]=1;
            lone[o]=rone[o]=mone[o]=one[o];
            lzr[o]=rzr[o]=mzr[o]=len[o]-one[o];
        }else{
            int m=l+r>>1;
            build(ls,l,m);
            build(rs,m+1,r);
            len[o]=len[ls]+len[rs];
            pu(o);
        }
    }
    void putone(int o,int l,int r){
        if(l>R||r<L)return;
        if(l>=L&&r<=R){
            To1(o);
        }else{
            pd(o);
            int m=l+r>>1;
            putone(ls,l,m);
            putone(rs,m+1,r);
            pu(o);
        }
    }
    void putzero(int o,int l,int r){
        if(l>R||r<L)return;
        if(l>=L&&r<=R){
            To0(o);
        }else{
            pd(o);
            int m=l+r>>1;
            putzero(ls,l,m);
            putzero(rs,m+1,r);
            pu(o);
        }
    }
    void reverse(int o,int l,int r){
        if(l>R||r<L)return;
        if(l>=L&&r<=R){
            Rs(o);
        }else{
            pd(o);
            int m=l+r>>1;
            reverse(ls,l,m);
            reverse(rs,m+1,r);
            pu(o);
        }
    }
    int getone(int o,int l,int r){
        if(l>R||r<L)return 0;
        if(l>=L&&r<=R){
            return one[o];
        }else{
            pd(o);
            int m=l+r>>1;
            return getone(ls,l,m)+getone(rs,m+1,r);
        }
    }
    int getMCO(int o,int l,int r){
        if(l>R||r<L)return 0;
        if(l>=L&&r<=R){
            return mone[o];
        }else{
            pd(o);
            int m=l+r>>1;
            int ml=getMCO(ls,l,m);
            int mr=getMCO(rs,m+1,r);
            int res=max(ml,mr);
            if(r<=R){
                res=max(res,lone[rs]+min(rone[ls],m+1-L));
            }else if(l>=L){
                res=max(res,rone[ls]+min(lone[rs],R-m));
            }else{
                res=max(res,min(rone[ls],m+1-L)+min(lone[rs],R-m));
            }
            return res;
        }
    }
}seg;


int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
    }
    seg.build(1,1,n);
    int op;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d",&op,&L,&R);
        L++;R++;
        if(op==0)seg.putzero(1,1,n);//置0
        else if(op==1)seg.putone(1,1,n);//置1
        else if(op==2)seg.reverse(1,1,n);//翻转
        else if(op==3)printf("%d\n",seg.getone(1,1,n));//区间和
        else printf("%d\n",seg.getMCO(1,1,n));//区间最大连续段
    }
}
```



### 线段树 - 带修线性基

```c++
//区间异或v 区间查询任意个数最大异或值
//用差分数组建立线段树
#include<bits/stdc++.h>
#define ls o<<1
#define rs o<<1|1
using namespace std;
const int M=50010;
const int N=31;
int n,m,a[M],b[M];

struct base{
    int a[N];
    base(){memset(a,0,sizeof a);}
    void add(int x){
        for(int i=N-1;i>=0&&x;--i)
            if(x>>i){
                if(!a[i]){
                    a[i]=x;
                    return;
                }
                x^=a[i];
            }
    }
    int get(int x){
        for(int i=N-1;i>=0;--i)x=max(x,x^a[i]);
        return x;
    }
    base merge(base t){
        for(int i=0;i<N;++i)t.add(a[i]);
        return t;
    }
}bs[M<<2];
int sum[M<<2];
int p,val,L,R;
void pu(int o){
    sum[o]=sum[ls]^sum[rs];
    bs[o]=bs[ls].merge(bs[rs]);
}
void build(int o,int l,int r){
    if(l==r){
        sum[o]=b[l];
        bs[o].add(sum[o]);
    }else{
        int m=l+r>>1;
        build(ls,l,m);
        build(rs,m+1,r);
        pu(o);
    }
}
int getSum(int o,int l,int r,int L,int R){
    if(l>R||r<L)return 0;
    if(l>=L&&r<=R)return sum[o];
    int m=l+r>>1;
    return getSum(ls,l,m,L,R)^getSum(rs,m+1,r,L,R);
}
void _add(int o,int l,int r){
    if(l==r){
        sum[o]^=val;
        bs[o]=base();
        bs[o].add(sum[o]);
    }else{
        int m=l+r>>1;
        if(p<=m)_add(ls,l,m);
        else _add(rs,m+1,r);
        pu(o);
    }
}
void add(int l,int r,int v){ //[l,r]区间异或v
    val=v;p=l;
    _add(1,1,n);
    p=r+1;
    if(p<=n)_add(1,1,n);
}
base _get(int o,int l,int r){
    if(l>R||r<L)return base();
    if(l>=L&&r<=R)return bs[o];
    int m=l+r>>1;
    return _get(ls,l,m).merge(_get(rs,m+1,r));
}
int get(int l,int r,int v){ //[l,r]内任意个数与v异或的最大值
    L=l+1;R=r;
    base t;
    if(L<=R)t=_get(1,1,n);
    t.add(getSum(1,1,n,1,l));
    return t.get(v);
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        b[i]=a[i-1]^a[i];//用差分数组建立线段树。
    }
    build(1,1,n);
    int op,l,r,v;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d %d",&op,&l,&r,&v);
        if(op==1) add(l,r,v);
        else printf("%d\n",get(l,r,v));
    }
}
```



### 可持久化数组

```c++
#include<bits/stdc++.h>
using namespace std;
const int M=1000010;
int n,m,cnt,a[M],rt[M];

struct nd{
    int l,r,x;
}tr[M<<4];

void build(int& o,int l,int r){
    o=++cnt;
    if(l==r)tr[o].x=a[l];
    else {
        int m=l+r>>1;
        build(tr[o].l,l,m);
        build(tr[o].r,m+1,r);
    }
}
int clone(int o){
    tr[++cnt]=tr[o];
    return cnt;
}
int p,val,ver;
void add(int& o,int l,int r){
    o=clone(o);
    if(l==r){
        tr[o].x=val;
    }else{
        int m=l+r>>1;
        if(p<=m)add(tr[o].l,l,m);
        else add(tr[o].r,m+1,r);
    }
}
int get(int o,int l,int r){
    if(l==r){
        return tr[o].x;
    }else{
        int m=l+r>>1;
        if(p<=m)return get(tr[o].l,l,m);
        else return get(tr[o].r,m+1,r);
    }
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    build(rt[0],1,n);
    int op;
    for(int i=1;i<=m;++i){
        scanf("%d %d",&ver,&op);
        rt[i]=rt[ver];
        if(op==1){
            scanf("%d %d",&p,&val);//修改
            add(rt[i],1,n);
        }else{
            scanf("%d",&p);
            printf("%d\n",get(rt[i],1,n));// 查询
        }
    }
}
```



### 可持久化线段树

```c++
//区间第k小
#include<bits/stdc++.h>
using namespace std;
const int M=200010;
int n,m,cnt,rt[M],a[M],rk[M],rnt;

struct nd{
    int x,l,r;
}tr[M<<5];

void pu(int o){
    tr[o].x=tr[tr[o].l].x+tr[tr[o].r].x;
}
int clone(int o){
    tr[++cnt]=tr[o];
    return cnt;
}
int p,val=1,k;
void add(int& o,int l,int r){
    o=clone(o);
    if(l==r){
        tr[o].x+=val;
    }else{
        int m=l+r>>1;
        if(p<=m)add(tr[o].l,l,m);
        else add(tr[o].r,m+1,r);
        pu(o);
    }
}
int get(int lo,int ro,int l,int r){
    if(l==r){
        return l;
    }else{
        int m=l+r>>1;
        if(tr[ro].x-tr[lo].x<k)return -1;
        int t=tr[tr[ro].l].x-tr[tr[lo].l].x;
        if(t>=k)return get(tr[lo].l,tr[ro].l,l,m);
        k-=t;
        return get(tr[lo].r,tr[ro].r,m+1,r);
    }
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i),rk[++rnt]=a[i];
    sort(rk+1,rk+1+n);
    rnt=unique(rk+1,rk+1+rnt)-rk-1;
    for(int i=1;i<=n;++i){
        a[i]=lower_bound(rk+1,rk+1+rnt,a[i])-rk;
    }

    for(int i=1;i<=n;++i){
        p=a[i];
        rt[i]=rt[i-1];
        add(rt[i],1,rnt);//建树 - 权值线段树
    }
    int l,r;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d",&l,&r,&k);
        printf("%d\n",rk[get(rt[l-1],rt[r],1,rnt)]);
    }
}
```

```c++
//求任意区间最小的不能组合出来的数。(p4587 主席树好题)
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
const int maxn=1e9;
int n,a[M],m;
int rt[M],cnt;

struct nd{
    int s,l,r;
}tr[M*50];

int clone(int o){
    tr[++cnt]=tr[o];
    return cnt;
}
int p,val;
void add(int& o,int l,int r){
    o=clone(o);
    if(l==r){
        tr[o].s+=p;
    }else{
        int m=l+r>>1;
        if(p<=m)add(tr[o].l,l,m);
        else add(tr[o].r,m+1,r);
        tr[o].s=tr[tr[o].l].s+tr[tr[o].r].s;
    }
}
int get(int lo,int ro,int l,int r){
    if(l>p||lo==0&&ro==0)return 0;
    if(r<=p)return tr[ro].s-tr[lo].s;
    int m=l+r>>1;
    return get(tr[lo].l,tr[ro].l,l,m)+get(tr[lo].r,tr[ro].r,m+1,r);
}

int main(){
    scanf("%d",&n);
    for(int i=1;i<=n;++i){
        scanf("%d",&p);
        rt[i]=rt[i-1];
        add(rt[i],1,maxn);
    }
    scanf("%d",&m);
    int x,y;
    for(int i=1;i<=m;++i){
        scanf("%d %d",&x,&y);
        int ans=0;p=0;
        while(ans>=p){
            p=ans+1;
            ans=get(rt[x-1],rt[y],1,maxn);
        }
        printf("%d\n",ans+1);
    }
}
```



### 可持久化并查集

```c++
#include<bits/stdc++.h>
using namespace std;
const int M=150010;
int n,m;

int rt[M<<1];
struct seg{
    int x[M<<4],ls[M<<4],rs[M<<4],cnt;
    int deep[M<<4];
    int dp;
    void build(int& o,int l,int r){
        o=++cnt;
        if(l==r){
            x[o]=l;
            deep[o]=1;
        }else{
            int m=l+r>>1;
            build(ls[o],l,m);
            build(rs[o],m+1,r);
        }
    }
    void init(){
        cnt=0;
        build(rt[cnt++],1,n);
    }
    int clone(int o){
        ++cnt;
        x[cnt]=x[o];
        ls[cnt]=ls[o];
        rs[cnt]=rs[o];
        return cnt;
    }
    void add(int& o,int l,int r,int& p,int& val){
        o=clone(o);
        if(l==r){
            x[o]=val;
        }else{
            int m=l+r>>1;
            if(p<=m)add(ls[o],l,m,p,val);
            else add(rs[o],m+1,r,p,val);
        }
    }
    void addDeep(int o,int l,int r,int p,int val){
        if(l==r){
            deep[o]=val;
        }else{
            int m=l+r>>1;
            if(p<=m)addDeep(ls[o],l,m,p,val);
            else addDeep(rs[o],m+1,r,p,val);
        }
    }
    int get(int o,int l,int r,int& p){
        if(l==r){
            dp=deep[o];
            return x[o];
        }
        int m=l+r>>1;
        if(p<=m)return get(ls[o],l,m,p);
        else return get(rs[o],m+1,r,p);
    }
    int Root(int ver,int p){//查询
        int rp=get(rt[ver],1,n,p);
        if(p==rp)return p;
        else return Root(ver,rp);
    }
    void OrderMerge(int ver,int p1,int p2){//按秩合并
        int dp1,dp2,r1,r2;
        r1=Root(ver,p1);dp1=dp;
        r2=Root(ver,p2);dp2=dp;
        if(r1==r2)return;
        if(dp1>dp2)add(rt[ver],1,n,r2,r1);
        else{
            add(rt[ver],1,n,r1,r2);
            if(dp1==dp2)addDeep(rt[ver],1,n,r2,dp2+1);
        }
    }
}S;

int main(){
    scanf("%d %d",&n,&m);
    S.init();
    int op,x1,x2;
    for(int i=1;i<=m;++i){
        scanf("%d",&op);
        if(op==1){
            scanf("%d %d",&x1,&x2);//合并
            rt[i]=rt[i-1];
            S.OrderMerge(i,x1,x2);
        }else if(op==2){
            scanf("%d",&x1);//跳转版本
            rt[i]=rt[x1];
        }else {
            scanf("%d %d",&x1,&x2);//查询
            rt[i]=rt[i-1];
            int r1=S.Root(i,x1);
            int r2=S.Root(i,x2);
            printf("%d\n",r1==r2);
        }
    }
}
```



### 可持久化字典树

```c++
//求最大异或值
#include<bits/stdc++.h>
using namespace std;
const int M=1000010;
const int N=25;
int n,m,a[M],sum[M];
int rt[M],tr[M<<4][2],num[M<<4],ed[M<<4],cnt;

int clone(int o){
    ++cnt;
    tr[cnt][0]=tr[o][0];
    tr[cnt][1]=tr[o][1];
    num[cnt]=num[o];
    return cnt;
}
void add(int& o,int x){//添加一个数x
    o=clone(o);
    int now=o;
    num[now]++;
    for(int i=N-1;i>=0;--i){
        int t=(x>>i)&1;
        tr[now][t]=clone(tr[now][t]);
        now=tr[now][t];
        num[now]++;
    }
    ed[now]=x;
}
int getMax(int lo,int ro,int ex){//lo-ro之间异或ex后的最大的数
    for(int i=N-1;i>=0;--i){
        int t=(ex>>i)&1;
        int x=num[tr[ro][t^1]]-num[tr[lo][t^1]];
        if(x>0){
            ro=tr[ro][t^1];
            lo=tr[lo][t^1];
        }else{
            ro=tr[ro][t^0];
            lo=tr[lo][t^0];
        }
    }
    return ed[ro];
}
ll getK(int o,ll ex,int k){//字典树中异或ex后k大数
    for(int i=N-1;i>=0;--i){
        int t=(ex>>i)&1;
        int x=num[tr[o][1^t]];
        if(x>=k)o=tr[o][1^t];
        else o=tr[o][0^t],k-=x;
    }
    return ed[o]^ex;
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        sum[i]=sum[i-1]^a[i];
        rt[i]=rt[i-1];
        add(rt[i],sum[i]);
    }
    char s[5];
    int x,l,r;
    for(int i=1;i<=m;++i){
        scanf("%s",s);
        if(s[0]=='A'){
            scanf("%d",&x);
            ++n;
            sum[n]=sum[n-1]^x;
            rt[n]=rt[n-1];
            add(rt[n],sum[n]);
        }else{
            scanf("%d %d %d",&l,&r,&x);
            printf("%d\n",getMax(rt[max(0,l-2)],rt[r-1],x^sum[n])^x^sum[n]);
        }
    }
}
```



### 线段树合并&分裂

```c++
#include<bits/stdc++.h>
const int M=200010;
using namespace std;
typedef long long ll;
ll n,m,a[M];

int rt[M],cnt;
ll x[M<<5];
int ls[M<<5],rs[M<<5];
void pd(int& o){
    if(!o)o=++cnt;
}
void pu(int& o){
    x[o]=x[ls[o]]+x[rs[o]];
}
void build(int& o,int l,int r){
    o=++cnt;
    if(l==r){
        x[o]=a[l];
    }else{
        int m=l+r>>1;
        build(ls[o],l,m);
        build(rs[o],m+1,r);
        pu(o);
    }
}
int L,R,p;
ll val,k;
void move(int& lo,int& ro,int l,int r){//分裂: lo->ro
    if(l>R||r<L)return;
    if(!lo)return;
    if(l>=L&&r<=R){
        swap(ro,lo);
    }else{
        pd(ro);
        int m=l+r>>1;
        move(ls[lo],ls[ro],l,m);
        move(rs[lo],rs[ro],m+1,r);
        pu(ro);
        pu(lo);
    }
}
void merge(int& lo,int& ro){//合并: lo->ro
    if(!lo)return;
    if(!ro)swap(lo,ro);
    else{
        x[ro]+=x[lo];
        x[lo]=0;
        merge(ls[lo],ls[ro]);
        merge(rs[lo],rs[ro]);
    }
}
void add(int& o,int l,int r){
    pd(o);
    if(l==r){
        x[o]+=val;
    }else{
        int m=l+r>>1;
        if(p<=m)add(ls[o],l,m);
        else add(rs[o],m+1,r);
        pu(o);
    }
}
ll get(int& o,int l,int r){
    if(l>R||r<L)return 0;
    pd(o);
    if(l>=L&&r<=R){
        return x[o];
    }else{
        int m=l+r>>1;
        return get(ls[o],l,m)+get(rs[o],m+1,r);
    }
}
int getKth(int& o,int l,int r,ll k){
    pd(o);
    if(l==r){
        return l;
    }else{
        if(x[o]<k)return -1;
        int m=l+r>>1;
        if(k<=x[ls[o]])return getKth(ls[o],l,m,k);
        else return getKth(rs[o],m+1,r,k-x[ls[o]]);
    }
}

int main(){
    scanf("%lld %lld",&n,&m);
    for(int i=1;i<=n;++i)scanf("%lld",a+i);
    int op,x,y;
    int top=1;
    build(rt[top],1,n);
    for(int i=1;i<=m;++i){
        scanf("%d",&op);
        if(op==0){
            scanf("%d %d %d",&x,&L,&R);
            move(rt[x],rt[++top],1,n);//分裂
        }else if(op==1){
            scanf("%d %d",&x,&y);
            merge(rt[y],rt[x]);//合并
        }else if(op==2){
            scanf("%d %lld %lld",&x,&val,&p);
            add(rt[x],1,n);//单点修改
        }else if(op==3){
            scanf("%d %d %d",&x,&L,&R);
            printf("%lld\n",get(rt[x],1,n));//区间查询
        }else {
            scanf("%d %lld",&x,&k);
            printf("%d\n",getKth(rt[x],1,n,k));//第k小
        }
    }
}
```



### 树套树 - 树状数组维护线段树

```c++
#include<bits/stdc++.h>
using namespace std;
const int M=600010;
int n,m,a[M],rt[M],rk[M],top;
//树状数组维护动态开点线段树 - 动态区间第k小
int x[M<<5],ls[M<<5],rs[M<<5],cnt;
int lowbit(int x){return x&-x;}
int p,L,R,val;
void pu(int& o){x[o]=x[ls[o]]+x[rs[o]];}
void pd(int& o){if(!o)o=++cnt;}
void _add(int& o,int l,int r){
    pd(o);
    if(l==r)x[o]+=val;
    else{
        int m=l+r>>1;
        if(p<=m)_add(ls[o],l,m);
        else _add(rs[o],m+1,r);
        pu(o);
    }
}
void add(int o,int p0,int val0){//修改权值线段树
    p=p0;val=val0;
    for(;o<=top;o+=lowbit(o))_add(rt[o],1,top);
}
int tmp[50][2];
int _getKth(int ln,int rn,int l,int r,int k){
    if(l==r)return rk[l];
    int lx=0,rx=0;
    for(int i=1;i<=ln;++i)lx+=x[ls[tmp[i][0]]];
    for(int i=1;i<=rn;++i)rx+=x[ls[tmp[i][1]]];
    int m=l+r>>1;
    if((rx-lx)>=k){
        for(int i=1;i<=ln;++i)tmp[i][0]=ls[tmp[i][0]];
        for(int i=1;i<=rn;++i)tmp[i][1]=ls[tmp[i][1]];
        return _getKth(ln,rn,l,m,k);
    }else{
        for(int i=1;i<=ln;++i)tmp[i][0]=rs[tmp[i][0]];
        for(int i=1;i<=rn;++i)tmp[i][1]=rs[tmp[i][1]];
        return _getKth(ln,rn,m+1,r,k-(rx-lx));
    }
}
int getKth(int lo,int ro,int k){//lo-ro区间第k小
    int ln=0,rn=0;
    for(;lo>0;lo-=lowbit(lo))tmp[++ln][0]=rt[lo];
    for(;ro>0;ro-=lowbit(ro))tmp[++rn][1]=rt[ro];
    return _getKth(ln,rn,1,top,k);
}

struct nd{
    int op,l,r,k;
}b[M];

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        rk[++top]=a[i];
    }
    char s[5];
    for(int i=1;i<=m;++i){
        scanf("%s",s);
        if(s[0]=='Q'){
            b[i].op=0;
            scanf("%d %d %d",&b[i].l,&b[i].r,&b[i].k);
        }else{
            b[i].op=1;
            scanf("%d %d",&b[i].l,&b[i].r);
            rk[++top]=b[i].r;
        }
    }
    sort(rk+1,rk+1+top);
    top=unique(rk+1,rk+1+top)-rk-1;
    for(int i=1;i<=n;++i){
        a[i]=lower_bound(rk+1,rk+1+top,a[i])-rk;
        add(i,a[i],1);
    }

    for(int i=1;i<=m;++i){
        if(b[i].op==0){
            printf("%d\n",getKth(b[i].l-1,b[i].r,b[i].k));
        }else{
            int t1=a[b[i].l];
            int t2=lower_bound(rk+1,rk+1+top,b[i].r)-rk;
            a[b[i].l]=t2;
            add(b[i].l,t1,-1);
            add(b[i].l,t2,1);
        }
    }
}
```



### 二维偏序 

二维偏序问题：离线问题树状数组，在线问题主席树。

二分套二维偏序问题为在线二维偏序，用主席树二维偏序（或整体二分）可解。

```c++
//归并排序求逆序对
#include<bits/stdc++.h>
using namespace std;
const int M=500010;
int a[M],b[M],n;
typedef long long ll;
ll ans;

void Bmerge(int l,int m,int r){
    int x=l,y=m+1,p=l;
    while(p<=r){
        if(y>r||a[x]<=a[y]&&x<=m)b[p++]=a[x++];
        else if(x>m||a[x]>a[y])b[p++]=a[y++],
        ans+=m+1-x;//<-
    }
    while(l<=r)a[l]=b[l],l++;
}

void Bsort(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    Bsort(l,m);
    Bsort(m+1,r);
    Bmerge(l,m,r);
}

int main(){
    scanf("%d",&n);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    Bsort(1,n);
    printf("%lld",ans);
}
```

```c++
//在线二维偏序(p4602)
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int M=100010;
const int N=100000;
int n,m;
struct nd{
    int d,p,l;
    bool operator < (const nd& t)const {
        return d<t.d;
    }
}a[M];
struct nd2{
    int l,r;
    ll x,p;
}tr[M*50];

int rt[M],cnt,val,p;
int clone(int o){
    tr[++cnt]=tr[o];
    return cnt;
}
void pu(int o){
    tr[o].x=tr[tr[o].l].x+tr[tr[o].r].x;
    tr[o].p=tr[tr[o].l].p+tr[tr[o].r].p;
}
void add(int &o,int l,int r){
    o=clone(o);
    if(l==r){
        tr[o].x+=a[p].l;
        tr[o].p=tr[o].x*l;
    }else{
        int m=l+r>>1;
        if(a[p].p<=m)add(tr[o].l,l,m);
        else add(tr[o].r,m+1,r);
        pu(o);
    }
}
ll get(int o,int l,int r,ll p){
    if(l==r)return min(p/l,tr[o].x);
    int m=l+r>>1;
    ll t=tr[tr[o].l].p;
    if(t<=p)return tr[tr[o].l].x+get(tr[o].r,m+1,r,p-t);
    else return get(tr[o].l,l,m,p);
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)
        scanf("%d %d %d",&a[i].d,&a[i].p,&a[i].l);
    sort(a+1,a+n+1);
    for(p=n;p>0;--p){
        rt[p]=rt[p+1];
        add(rt[p],1,N);
    }

    ll x,y;
    for(int i=1;i<=m;++i){
        scanf("%lld %lld",&x,&y);
        int l=0,r=n;
        while(l<r){
            int mid=(l+r+1)>>1;
            if(get(rt[mid],1,N,x)>=y)l=mid;
            else r=mid-1;
        }
        if(l)printf("%d\n",a[l].d);
        else printf("-1\n");
    }
}
```



### 三维偏序

#### 单向开区间

```c++
//cdq分治模板题 单向开区间三维偏序
#include<bits/stdc++.h>
using namespace std;
const int M=200010;
int n,top,m,ans[M];

struct nd{
    int x,y,z,cnt,ans;
}a[M],b[M];
bool cmpy(const nd& t1,const nd& t2){
    if(t1.y==t2.y)return t1.z<t2.z;
    return t1.y<t2.y;
}
bool cmpx(const nd& t1,const nd& t2){
    if(t1.x==t2.x)return cmpy(t1,t2);
    return t1.x<t2.x;
}

struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int v){
        for(;x<=m;x+=lowbit(x))a[x]+=v;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))res+=a[x];
        return res;
    }
}bt;

void cdq(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    cdq(l,m);
    cdq(m+1,r);
    sort(a+l,a+m+1,cmpy);
    sort(a+m+1,a+r+1,cmpy);

    int x=l,y=m+1;
    while(y<=r){
        if(x<=m&&a[x].y<=a[y].y){
            bt.add(a[x].z,a[x].cnt);
            x++;
        }else{
            a[y].ans+=bt.get(a[y].z);
            y++;
        }
    }
    for(;l<x;++l)bt.add(a[l].z,-a[l].cnt);
}


int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d %d %d",&b[i].x,&b[i].y,&b[i].z);
    }
    sort(b+1,b+1+n,cmpx);
    int t=0;
    for(int i=1;i<=n;++i){
        t++;
        if(b[i].x!=b[i+1].x||b[i].y!=b[i+1].y||b[i].z!=b[i+1].z){
            a[++top].x=b[i].x;
            a[top].y=b[i].y;
            a[top].z=b[i].z;
            a[top].cnt=t;
            t=0;
        }
    }
	//三维偏序
    cdq(1,top);
    for(int i=1;i<=top;++i)
        ans[a[i].ans+a[i].cnt-1]+=a[i].cnt;
    for(int i=0;i<n;++i)printf("%d\n",ans[i]);
}
```

```c++
// 矩形区间三维偏序
#include<bits/stdc++.h>
using namespace std;
const int M=400010;
int n,k,top;
int rk[M];
long long ans;
struct nd{
    int x,r,q;
    int vl,vr,ql,qr;
    bool operator < (const nd& t) const {
        return r>t.r;
    }
}a[M];

bool cmp1(nd& t1,nd& t2){
    return t1.q<t2.q;
}

struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int v){
        for(;x<=top;x+=lowbit(x))
            a[x]+=v;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))
            res+=a[x];
        return res;
    }
}bt;

void cdq(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    cdq(l,m);
    cdq(m+1,r);
    sort(a+l,a+m+1,cmp1);
    sort(a+m+1,a+r+1,cmp1);

    //q的匹配范围为长度为2k的定长区间，所以可以直接分治q;
    //或者将x的不定长匹配区间拆开（左、中、右），分治x。
    //这里分治q
    int pl=l,pr=l,x=m+1;
    while(x<=r){
        while(pr<=m&&a[pr].q<=a[x].q+k)bt.add(a[pr].x, 1),pr++;
        while(pl<=m&&a[pl].q< a[x].q-k)bt.add(a[pl].x,-1),pl++;
        ans+=bt.get(a[x].vr)-bt.get(a[x].vl-1);
        x++;
    }
    while(pl<pr)bt.add(a[pl].x,-1),pl++;
}

int main(){
    scanf("%d %d",&n,&k);
    for(int i=1;i<=n;++i){
        scanf("%d %d %d",&a[i].x,&a[i].r,&a[i].q);
        rk[++top]=a[i].x;
        rk[++top]=a[i].x+a[i].r;
        rk[++top]=a[i].x-a[i].r;
    }
    sort(rk+1,rk+1+top);
    top=unique(rk+1,rk+1+top)-rk-1;
    for(int i=1;i<=n;++i){
        a[i].vl=lower_bound(rk+1,rk+top+1,a[i].x-a[i].r)-rk;
        a[i].vr=lower_bound(rk+1,rk+top+1,a[i].x+a[i].r)-rk;
        a[i].x=lower_bound(rk+1,rk+top+1,a[i].x)-rk;
    }
    sort(a+1,a+1+n);
    cdq(1,n);
    printf("%lld",ans);
}
```

#### 三维数点

> 问题模型：
>
> 已知：M=1e5，a[M], b[M], c[M];
> M次询问：
> 每次给出: l1, r1, l2, r2, l3, r3;
> 每次询问: 满足 l1<a[x]<r1, l2<b[x]<r2, l3<c[x]<r3 的x的**个数**。

------

> cdq分治：
>
> 将**数据**和**询问**放入同一个数组，按照**前两维**进行排序；每1个第一二维的区间询问拆解为4个单向询问。
> 每一层归并时，只计算左边的数据对右边的询问的贡献：
> 	1 整体排序，左边整体小于右边整体，直接筛选出第一维贡献。
> 	2 归并排序，左边和右边内部分别有序，双指针筛选第二维贡献。
> 	3 树状数组，此时第一二维已经被筛选掉，所以o(logN)维护&查询这一层的答案。
> NlogN次归并，累计的答案即为最终结果。
>
> **注意：先插入数据，再插入询问，排序函数使用 stable_sort ，并且只能针对前两维进行排序，以保证数据均出现在询问的左边。**

```c++
//https://www.luogu.com.cn/problem/P4396
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
int n,m,a[M],pre[M],lst[M],ans1[M],ans2[M];

struct bit{
    int a[M];
    int lowbit(int x){return x&-x;};
    void add(int x,int val){
        for(;x<M;x+=lowbit(x))
            a[x]+=val;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))
            res+=a[x];
        return res;
    }
}tr;

struct nd{
    int x,y,z,f,id;
    bool operator < (const nd& t)const{
        return y<t.y;
    }
}b[M*5];
int cnt;
bool cmp(const nd& t1,const nd& t2){
    if(t1.x==t2.x)return t1<t2;
    return t1.x<t2.x;
}

void cdq(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    cdq(l,m);
    cdq(m+1,r);

    int x=l,y=m+1,t=0;
    while(y<=r){
        if(x<=m&&b[x].y<=b[y].y){
            if(b[x].f==0)tr.add(b[x].z+1,1),t++;
            x++;
        }else{
            if(b[y].f!=0){
                ans1[b[y].id]+=t*b[y].f;
                ans2[b[y].id]+=tr.get(b[y].z+1)*b[y].f;
            }
            y++;
        }
    }
    for(int i=l;i<x;++i)
        if(b[i].f==0)tr.add(b[i].z+1,-1);
    stable_sort(b+l,b+r+1);
    //归并排序中stable_sort比sort快。
}

int main(){
    scanf("%d %d" ,&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        pre[i]=lst[a[i]];
        lst[a[i]]=i;
        b[++cnt]=(nd){i,a[i],pre[i],0,0};
    }
    int l,r,A,B;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d %d",&l,&r,&A,&B);//询问拆解
        b[++cnt]=(nd){r,B,l-1,1,i};
        b[++cnt]=(nd){l-1,A-1,l-1,1,i};
        b[++cnt]=(nd){r,A-1,l-1,-1,i};
        b[++cnt]=(nd){l-1,B,l-1,-1,i};
    }
    stable_sort(b+1,b+cnt+1,cmp);
    cdq(1,cnt);
    for(int i=1;i<=m;++i)printf("%d %d\n",ans1[i],ans2[i]);
}
```



### 整体二分

整体二分与cdq分治有异曲同工之妙。

二分答案的过程中即为对所有询问分治的过程。

```c++
//https://www.luogu.com.cn/problem/P1527
/*
N=500,给出N*N的矩阵
1e5次询问，每次询问一个矩形子区间中第k大的值。
*/
#include<bits/stdc++.h>
using namespace std;
const int M=60010;
const int N=510;
int n,m,top,ans[M];

//要分治 先排序
struct nd{
    int x,y,v;
    bool operator < (const nd& t)const {
        return v<t.v;
    }
}a[N*N];
struct Q{
    int x1,y1,x2,y2,k;
    int id;
};
typedef vector<Q> VQ;

struct BIT2{
    int a[N][N];
    int lowbit(int x){return x&-x;}
    void add(int x,int y,int v){
        for(;x<=n;x+=lowbit(x))
            for(int t=y;t<=n;t+=lowbit(t))
                a[x][t]+=v;
    }
    int sum(int x,int y){
        int res=0;
        for(;x>0;x-=lowbit(x))
            for(int t=y;t>0;t-=lowbit(t))
                res+=a[x][t];
        return res;
    }
    int get(int x1,int y1,int x2,int y2){
        int res=sum(x1-1,y1-1)+sum(x2,y2);
        res-=sum(x2,y1-1)+sum(x1-1,y2);
        return res;
    }
}tr;

void cdq(int l,int r,VQ vc){
    if(vc.size()==0)return;
    if(l==r){
        for(auto v:vc)ans[v.id]=a[l].v;
        return;
    }
    int m=l+r>>1;
    for(int i=l;i<=m;++i)
        tr.add(a[i].x,a[i].y,1);
    VQ v1,v2;
    for(auto v:vc){
        int t=tr.get(v.x1,v.y1,v.x2,v.y2);
        if(v.k<=t)v1.push_back(v);
        else{
            v.k-=t;
            v2.push_back(v);
        }
    }
    for(int i=l;i<=m;++i)
        tr.add(a[i].x,a[i].y,-1);

    cdq(l,m,v1);
    cdq(m+1,r,v2);
}

int main(){
    scanf("%d %d" ,&n,&m);
    for(int i=1;i<=n;++i){
        for(int j=1;j<=n;++j){
            a[++top].x=i;
            a[top].y=j;
            scanf("%d",&a[top].v);
        }
    }
    sort(a+1,a+1+top);
    VQ vc;
    Q t;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d %d %d",&t.x1,&t.y1,&t.x2,&t.y2,&t.k);
        t.id=i;
        vc.push_back(t);
    }
    //cdq离线处理
    cdq(1,top,vc);
    for(int i=1;i<=m;++i){
        printf("%d\n",ans[i]);
    }
}
```



### Splay 

```c++
//单点插入、删除 前驱 后继 第k大 求排名
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
const int inf=1<<30;

struct Splay{
    int x[M]={-inf},fa[M],ls[M],rs[M],sum[M],num[M];
    int& root=rs[0];
    int cnt;
    void pu(int o){
        sum[o]=sum[ls[o]]+sum[rs[o]]+num[o];
    }
    int getdir(int o){
        return ls[fa[o]]==o?0:1;
    }
    void connect(int o,int f,int dir){
        fa[o]=f;
        if(dir)rs[f]=o;
        else ls[f]=o;
    }
    void rotate(int o){
        int f=fa[o];
        int ff=fa[f];
        int dirf=getdir(f);
        int dirx=getdir(o);
        int t=(dirx^1)?rs[o]:ls[o];
        connect(t,f,dirx);connect(f,o,dirx^1);connect(o,ff,dirf);
        pu(f);pu(o);
    }
    void splay(int o,int to){
        to=fa[to];
        while(fa[o]!=to){
            int f=fa[o];
            if(fa[f]==to)rotate(o);
            else if(getdir(o)==getdir(f)){
                rotate(f);
                rotate(o);
            }else{
                rotate(o);
                rotate(o);
            }
        }
    }
    int newpoint(int v,int f){
        x[++cnt]=v;
        fa[cnt]=f;
        sum[cnt]=num[cnt]=1;
        return cnt;
    }
    int find(int v){
        int now=root;
        while(1){
            if(x[now]==v){
                splay(now,root);
                return now;
            }
            now=v<x[now]?ls[now]:rs[now];
            if(now==0)return 0;
        }
    }
    ///---
    void add(int v){
        int now=root;
        while(1){
            if(now)sum[now]++;
            if(v==x[now]){
                num[now]++;
                break;
            }
            int& t=v<x[now]?ls[now]:rs[now];
            if(t==0){
                now=t=newpoint(v,now);
                break;
            }
            now=t;
        }
        splay(now,root);
    }
    void remove(int v){
        int o=find(v);
        if(o==0)return;
        if(num[o]>1){
            num[o]--;
            sum[o]--;
            return;
        }
        if(ls[o]==0){
            root=rs[o];
            fa[root]=0;
        }else{
            int l=ls[o];
            while(rs[l])l=rs[l];
            splay(l,ls[o]);
            connect(rs[o],l,1);
            connect(l,0,1);
            pu(l);
        }
    }
    int getRank(int v){
        int res=0,o=root;
        while(1){
            if(o==0)break;
            if(x[o]==v){
                res+=sum[ls[o]]+1;
                break;
            }
            if(v<x[o])o=ls[o];
            else {
                res+=sum[ls[o]]+num[o];
                o=rs[o];
            }
        }
        if(o)splay(o,root);
        return res;
    }
    int rank(int p){
        if(p<=0||p>sum[root])return -inf;
        int o=root;
        while(1){
            int m=sum[o]-sum[rs[o]];
            if(p>sum[ls[o]] && p<=m)break;
            if(p<m)o=ls[o];
            else {
                p-=m;
                o=rs[o];
            }
        }
        splay(o,root);
        return x[o];
    }
    int upper(int v){
        int o=root,res=inf;
        while(o){
            if(x[o]>v && x[o]<res)res=x[o];
            if(v<x[o])o=ls[o];
            else o=rs[o];
        }
        return res;
    }
    int lower(int v){
        int o=root,res=-inf;
        while(o){
            if(x[o]<v && x[o]>res)res=x[o];
            if(v>x[o])o=rs[o];
            else o=ls[o];
        }
        return res;
    }
};
```



### 莫队 - 静态区间众数

```c++
//值域1e9,离线求任意区间最多的数的个数
#include<bits/stdc++.h>
using namespace std;
const int M=500010;
const int maxn=1e9+7;
int n,m,a[M],r[M];
struct nd{
    int l,r,L,id;
    bool operator < (const nd& t)const {
        if(L!=t.L)return L<t.L;
        return (r<t.r);
    }
}b[M];

int now,cnt[M],ans[M],num[M];
void add(int t){
    num[cnt[t]]--;
    cnt[t]++;
    num[cnt[t]]++;
    if(cnt[t]>now)now=cnt[t];
}
void del(int t){
    num[cnt[t]]--;
    cnt[t]--;
    num[cnt[t]]++;
    if(cnt[t]+1==now)
        if(num[now]==0)now--;
}

int main(){
    scanf("%d %d" ,&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        r[i]=a[i];
    }
    sort(r+1,r+1+n);
    int nr=unique(r+1,r+1+n)-r-1;
    for(int i=1;i<=n;++i)a[i]=lower_bound(r+1,r+1+nr,a[i])-r;

    int s=sqrt(n);
    for(int i=1;i<=m;++i){
        scanf("%d %d",&b[i].l,&b[i].r);
        b[i].id=i;
        b[i].L=b[i].l/s;
    }
    sort(b+1,b+1+m);
    int l=1,r=0;
    for(int i=1;i<=m;++i){
        int ll=b[i].l,rr=b[i].r;
        while(r<rr)add(a[++r]);
        while(l>ll)add(a[--l]);
        while(r>rr)del(a[r--]);
        while(l<ll)del(a[l++]);
        ans[b[i].id]=now;
    }
    for(int i=1;i<=m;i++)printf("%d\n",ans[i]);
}
```



### 不删除莫队

```c++
//https://www.luogu.com.cn/problem/AT1219
//求区间内 x*cnt[x] 的最大值。
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
typedef long long ll;
int n,m,a[M],rk[M];

struct nd{
    int l,r,L,id;
    bool operator <(const nd& t)const{
        if(L!=t.L)return L<t.L;
        return r<t.r;
    }
}b[M],c[M];
int cntb,cntc;

int cnt[M],cnt0[M];
ll now,now0,ans[M];

void force(int l,int r){//块内暴力
    now=0;
    for(int i=l;i<=r;++i)cnt[a[i]]=0;
    for(int i=l;i<=r;++i){
        cnt[a[i]]++;
        now=max(now,1ll*cnt[a[i]]*rk[a[i]]);
    }
}
void add(int x){
    cnt[x]++;
    now=max(now,1ll*cnt[x]*rk[x]);
}
void addtmp(int x){ // 临时修改
    cnt0[x]++;
    now0=max(now0,1ll*(cnt[x]+cnt0[x])*rk[x]);
}
void clear(int l,int r){ // 回收临时修改
    now0=0;
    for(int i=l;i<=r;++i)cnt0[a[i]]=0;
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i),rk[i]=a[i];
    sort(rk+1,rk+1+n);
    for(int i=1;i<=n;++i)
        a[i]=lower_bound(rk+1,rk+1+n,a[i])-rk;

    int l,r,s=n/sqrt(m);
    for(int i=1;i<=m;++i){
        scanf("%d %d",&l,&r);
        if(l/s==r/s)c[++cntc]=(nd){l,r,0,i};
        else b[++cntb]=(nd){l,r,l/s+1,i};
    }
    for(int i=1;i<=cntc;++i){
        force(c[i].l,c[i].r);
        ans[c[i].id]=now;
    }
    sort(b+1,b+1+cntb);
    for(int i=1;i<=cntb;++i){
        l=b[i].L*s;
        if(b[i].L!=b[i-1].L){
            memset(cnt,0,sizeof cnt);
            now=0;
            r=l-1;
        }
        while(r<b[i].r)add(a[++r]);
        while(l>b[i].l)addtmp(a[--l]);
        ans[b[i].id]=max(now,now0);
        clear(l,b[i].L*s);
    }
    for(int i=1;i<=m;++i)printf("%lld\n",ans[i]);
}
```



### 带修莫队

```c++
#include<bits/stdc++.h>
using namespace std;
const int N=200010;
const int M=1000010;
int n,m,a[N],ts[N],b[N],ans[N];
int top,cnt;

struct nd{
    int l,r,L,R,t,id;
    bool operator <(const nd& v)const{
        if(L!=v.L)return L<v.L;
        if(R!=v.R)return (R<v.R);
        return (t<v.t);
    }
}q[N];

int now,num[M];
void add(int x){
    if(num[x]==0)now++;
    num[x]++;
}
void del(int x){
    num[x]--;
    if(num[x]==0)now--;
}
void update(int x,int l,int r){
    int i=ts[x];
    if(i>=l&&i<=r){
        add(b[x]);
        del(a[i]);
    }
    swap(a[i],b[x]);
}

int main(){
    char c[10];
    int x,y;
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    for(int i=1;i<=m;++i){
        scanf("%s %d %d",c,&x,&y);
        if(c[0]=='R'){
            ts[++top]=x;
            b[top]=y;
        }else{
            ++cnt;
            q[cnt]=(nd){x,y,0,0,top,cnt};
        }
    }

    int s=pow(n,0.666);
    for(int i=1;i<=cnt;++i){
        q[i].L=q[i].l/s;
        q[i].R=q[i].r/s;
    }
    sort(q+1,q+1+cnt);
    int l=1,r=0,t=0;
    for(int i=1;i<=cnt;++i){
        int x=q[i].l,y=q[i].r,z=q[i].t;
        while(r<y)add(a[++r]);
        while(l>x)add(a[--l]);
        while(r>y)del(a[r--]);
        while(l<x)del(a[l++]);
        while(t<z)update(++t,l,r);
        while(t>z)update(t--,l,r);
        ans[q[i].id]=now;
    }
    for(int i=1;i<=cnt;++i)printf("%d\n",ans[i]);
}
```



### 树上莫队

```c++
//https://www.luogu.com.cn/problem/SP10707
#include<bits/stdc++.h>
using namespace std;
const int M=200010;
int n,m,a[M],rk[M],ans[M];
struct nd_q{
    int x,id;
};
vector<int> g[M];
vector<nd_q> g_q[M];
struct pt{
    int l,r,L,id,lca;
    bool operator < (const pt& t)const{
        if(L!=t.L)return L<t.L;
        return (r<t.r);
    }
}b[M];

int root[M],vis[M];
int find(int x){return root[x]==x?x:root[x]=find(root[x]);}
void merge(int x,int y){root[find(x)]=find(y);}
void tarjan(int o){
    vis[o]=1;
    for(auto v:g[o]){
        if(vis[v]==0){
            tarjan(v);
            merge(v,o);
        }
    }
    for(auto v:g_q[o]){
        if(vis[v.x]==1)b[v.id].lca=find(v.x);
    }
}//树上莫队会用到lca

int cnt,st[M],ed[M],pot[M];
void build(int o,int fa){
    st[o]=++cnt;pot[cnt]=o;
    for(auto v:g[o]){
        if(fa!=v)build(v,o);
    }
    ed[o]=++cnt;pot[cnt]=o;
}//欧拉序展开

int now,use[M],num[M];
void add(int x){
    use[x]^=1;
    if(use[x]){
        if(++num[a[x]] == 1)now++;
    }else{
        if(--num[a[x]] == 0)now--;
    }
}

int main(){
    for(int i=1;i<M;++i)root[i]=i;
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        rk[i]=a[i];
    }
    sort(rk+1,rk+1+n);
    for(int i=1;i<=n;++i)a[i]=lower_bound(rk+1,rk+1+n,a[i])-rk;
    int x,y;
    for(int i=1;i<n;++i){
        scanf("%d %d",&x,&y);
        g[x].push_back(y);
        g[y].push_back(x);
    }
    for(int i=1;i<=m;++i){
        scanf("%d %d",&b[i].l,&b[i].r);
        b[i].id=i;
        g_q[b[i].l].push_back({b[i].r,i});
        g_q[b[i].r].push_back({b[i].l,i});
    }
    build(1,0);
    tarjan(1);
    int s=sqrt(n*2);
    for(int i=1;i<=m;++i){
        int& x=b[i].l;
        int& y=b[i].r;
        if(st[x]>st[y])swap(x,y);
        if(b[i].lca==x){
            b[i].lca=0;
            b[i].l=st[x];
            b[i].r=st[y];
        }else{
            b[i].l=ed[x];
            b[i].r=st[y];
        }
        b[i].L=b[i].l/s+1;
    }
    sort(b+1,b+1+m);
    int l=1,r=0;
    for(int i=1;i<=m;++i){
        while(r<b[i].r)add(pot[++r]);
        while(l>b[i].l)add(pot[--l]);
        while(r>b[i].r)add(pot[r--]);
        while(l<b[i].l)add(pot[l++]);
        if(b[i].lca)add(b[i].lca);
        ans[b[i].id]=now;
        if(b[i].lca)add(b[i].lca);
    }
    for(int i=1;i<=m;++i)printf("%d\n",ans[i]);
}
```



### 莫队二次离线

```c++
//https://www.luogu.com.cn/problem/P4887
//每次询问满足 ai^aj 二进制下有 K 个1的二元组个数。
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int M=100010;
int n,m,k,s,a[M];
struct nd{
    int l,r,L,id;
    ll ans_delta;
    bool operator < (const nd& t)const{
        if(L!=t.L)return L<t.L;
        return r<t.r;
    }
}b[M];
struct nd_q{
    int l,r,id,f;
};
vector<nd_q> q1[M];

ll ans[M];
int k1[M],top,pre_eff[M],cnt[M];
void pre_init(){
    for(int i=0;i<16384;++i){
        if(__builtin_popcount(i)==k)k1[++top]=i;
    }
    for(int i=1;i<=n;++i){
        pre_eff[i]=cnt[a[i]];
        for(int j=1;j<=top;++j)cnt[k1[j]^a[i]]++;
    }
}

int main(){
    scanf("%d %d %d",&n,&m,&k);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    s=sqrt(n);
    for(int i=1;i<=m;++i){
        scanf("%d %d",&b[i].l,&b[i].r);
        b[i].id=i;
        b[i].L=b[i].l/s+1;
    }
    if(k>14){
        for(int i=1;i<=m;++i)puts("0");
        return 0;
    }
    pre_init();
    sort(b+1,b+1+m);
    //莫队，拆解答案 或 计算部分答案。
    int l=1,r=0;
    for(int i=1;i<=m;++i){
        if(l<b[i].l)q1[r].push_back({l,b[i].l-1,i,-1});
        while(l<b[i].l)b[i].ans_delta+=pre_eff[l++];
        if(l>b[i].l)q1[r].push_back({b[i].l,l-1,i,1});
        while(l>b[i].l)b[i].ans_delta-=pre_eff[--l];
        if(r<b[i].r)q1[l-1].push_back({r+1,b[i].r,i,-1});
        while(r<b[i].r)b[i].ans_delta+=pre_eff[++r];
        if(r>b[i].r)q1[l-1].push_back({b[i].r+1,r,i,1});
        while(r>b[i].r)b[i].ans_delta-=pre_eff[r--];
    }
    // 再次离线计算答案
    memset(cnt,0,sizeof cnt);
    for(int i=1;i<=n;++i){
        for(int j=1;j<=top;++j)cnt[k1[j]^a[i]]++;
        for(auto v:q1[i])
            for(int j=v.l;j<=v.r;++j){
                if(k==0&&j<=i)b[v.id].ans_delta+=v.f*(cnt[a[j]]-1);
                else b[v.id].ans_delta+=v.f*cnt[a[j]];
            }
    }
    //差分还原
    for(int i=1;i<=m;++i)b[i].ans_delta+=b[i-1].ans_delta;
    for(int i=1;i<=m;++i)ans[b[i].id]=b[i].ans_delta;
    for(int i=1;i<=m;++i)printf("%lld\n",ans[i]);
    return 0;
}
```

```c++
//https://www.luogu.com.cn/problem/P5047
//询问任意区间逆序对
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int M=100010;
int n,m,a[M],s,rk[M];

struct nd{
    int l,r,L,id;
    ll ans;
    bool operator < (const nd& t)const{
        if(L!=t.L)return L<t.L;
        if(L&1)return r<t.r;
        else return r>t.r;
    }
}b[M];
struct Q{int l,r,f,id;};
vector<Q> qs[M];

struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int val){
        for(;x<=n;x+=lowbit(x))a[x]+=val;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))res+=a[x];
        return res;
    }
}tr;

struct BT{// o(1)查询
    int k[M],as[M],aq[M],s,ns;
    void init(){//!!!
        s=sqrt(n);ns=n/s+1;
        for(int i=1;i<=n;++i)k[i]=i/s+1;
    }
    void add(int x,int val){
        for(int i=k[x]+1;i<=ns;++i)as[i]+=val;
        for(int x0=x;x<=n&&k[x]==k[x0];++x)aq[x]+=val;
    }
    int get(int x){
        return as[k[x]]+aq[x];
    }
}bt;

int pre[M],suf[M];
ll ans[M];

int main(){
    scanf("%d %d",&n,&m);
    s=sqrt(n);
    for(int i=1;i<=n;++i)scanf("%d",a+i),rk[i]=a[i];
    sort(rk+1,rk+1+n);
    for(int i=1;i<=n;++i)
        a[i]=lower_bound(rk+1,rk+1+n,a[i])-rk;
    for(int i=1;i<=m;++i){
        scanf("%d %d",&b[i].l,&b[i].r);
        b[i].L=b[i].l/s+1;
        b[i].id=i;
    }
    sort(b+1,b+1+m);
    for(int i=1;i<=n;++i){
        pre[i]=tr.get(n)-tr.get(a[i]);
        tr.add(a[i],1);
    }
    memset(tr.a,0,sizeof tr.a);
    for(int i=n;i>0;--i){
        suf[i]=tr.get(a[i]-1);
        tr.add(a[i],1);
    }
	//将转移分为两部分：固定部分和移动部分，分开处理（if和while）。
    int l=1,r=0;
    for(int i=1;i<=m;++i){
        int ll=b[i].l,rr=b[i].r;
        if(l<ll)qs[r+1].push_back({l,ll-1,1,i});
        while(l<ll)b[i].ans-=suf[l++];
        if(l>ll)qs[r+1].push_back({ll,l-1,-1,i});
        while(l>ll)b[i].ans+=suf[--l];
        if(r<rr)qs[l-1].push_back({r+1,rr,-1,i*-1});
        while(r<rr)b[i].ans+=pre[++r];
        if(r>rr)qs[l-1].push_back({rr+1,r,1,i*-1});
        while(r>rr)b[i].ans-=pre[r--];
    }
	//对于逆序对问题，分成向左向右两部分，分开处理。
    bt.init();
    for(int i=1;i<=n;++i){
        bt.add(a[i],1);
        for(auto v:qs[i]){
            if(v.id<0)
            for(int j=v.l;j<=v.r;++j)
                b[v.id*-1].ans+=v.f*(bt.get(n)-bt.get(a[j]));
        }
    }
    memset(bt.aq,0,sizeof bt.aq);
    memset(bt.as,0,sizeof bt.as);
    for(int i=n;i>0;--i){
        bt.add(a[i],1);
        for(auto v:qs[i]){
            if(v.id>0)
            for(int j=v.l;j<=v.r;++j)b[v.id].ans+=v.f*bt.get(a[j]-1);
        }
    }
    for(int i=1;i<=m;++i)b[i].ans+=b[i-1].ans;
    for(int i=1;i<=m;++i)ans[b[i].id]=b[i].ans;
    for(int i=1;i<=m;++i)printf("%lld\n",ans[i]);
}
```





## 树上问题

### 直径

```c++
int dis1,dis2,id;

//两次dfs 
void dfs1(int o,int fa,int deep){
    if(deep>=dis1){
        id=o;
        dis1=deep;
    }
    for(auto v:vc[o]){
        if(v.x!=fa)dfs1(v.x,o,deep+1);
    }
}

//树上dp 可带负权边
void dfs2(int o,int fa){
    for(auto& v:vc[o]){
        if(v.x!=fa){
            dfs2(v.x,o);
            if(id==v.x){
                id=o;
                v.w=-1;
            }
        }
    }
}

```



### 重心

```c++
//删除重心后所得的所有子树，节点数不超过原树的1/2，一棵树最多有两个重心
int n;
vector<int> g[M];

//重心为rt
int sz[M],rt;
void dfs1(int o,int fa){
    sz[o]=1;
    int maxsz=0;
    for(auto v:g[o]){
        if(v!=fa){
            dfs1(v,o);
            sz[o]+=sz[v];
            maxsz=max(maxsz,sz[v]);
        }
    }
    maxsz=max(maxsz,n-sz[o]);
    if(maxsz<=n/2)rt=o;
}
```



### LCA

```c++
//离线lca
vector<int> vc[M];//树
int root[M],vis[M];
int find(int x){
    return root[x]==x?x:(root[x]=find(root[x]));
}
void merge(int x,int y){
    root[find(x)]=find(y);
}

struct nd{int v,id;};
vector<nd> q[M];//查询

void tarjan(int o){
    vis[o]=1;
    for(auto v:vc[o]){
        if(vis[v]==0){
            tarjan(v);
            merge(v,o);
        }
    }
    for(auto v:q[o]){
        if(vis[v.v]==1){
            ans[v.id]=find(v.v);
        }
    }
}

//在线倍增
"见红书  p74"
    
//树剖LCA
vector<int> g[M];
int deep[M],fa[M],son[M],sz[M];
int top[M];

void dfs1(int o,int f,int d){
    fa[o]=f;
    deep[o]=d;
    sz[o]=1;
    int maxsz=0;
    for(auto v:g[o]){
        if(v!=f){
            dfs1(v,o,d+1);
            sz[o]+=sz[v];
            if(sz[v]>maxsz){
                son[o]=v;
                maxsz=sz[v];
            }
        }
    }
}

void dfs2(int o,int t){
    top[o]=t;
    if(son[o]==0)return;
    dfs2(son[o],t);
    for(auto v:g[o]){
        if(v!=fa[o]&&v!=son[o]){
            dfs2(v,v);
        }
    }
}

int lca(int x,int y){
    while(top[x]!=top[y]){
        if(deep[top[x]]<deep[top[y]])swap(x,y);
        x=fa[top[x]];
    }
    if(deep[x]<deep[y])return x;
    else return y;
}
```



### 树链剖分

```c++
//p3384 树链剖分模板
#include<bits/stdc++.h>
#define ls o<<1
#define rs o<<1|1
using namespace std;
const int M=100010;
typedef long long ll;
ll n,m,rt,mod,a[M],b[M];
vector<int> g[M];

//tree
ll deep[M],fa[M],sz[M],son[M];
ll cnt,id[M],top[M];

void dfs1(int o,int f,int d){
    fa[o]=f;
    deep[o]=d;
    sz[o]=1;
    int maxson=0;
    for(auto v:g[o]){
        if(v!=f){
            dfs1(v,o,d+1);
            sz[o]+=sz[v];
            if(sz[v]>maxson){
                maxson=sz[v];
                son[o]=v;
            }
        }
    }
}
void dfs2(int o,int topf){
    id[o]=++cnt;
    b[cnt]=a[o];
    top[o]=topf;
    if(!son[o])return;
    dfs2(son[o],topf);
    for(auto v:g[o])
        if(id[v]==0)dfs2(v,v);
}

// seg
ll L,R,val;
ll x[M<<2],lz[M<<2],len[M<<2];

void do_1(int o,ll val){
    x[o]=(x[o]+val*len[o])%mod;
    lz[o]=(lz[o]+val)%mod;
}
void pu(int o){
    x[o]=(x[ls]+x[rs])%mod;
}
void pd(int o){
    if(lz[o]){
        do_1(ls,lz[o]);
        do_1(rs,lz[o]);
        lz[o]=0;
    }
}
void build(int o,int l,int r){
    if(l==r){
        x[o]=b[l];
        len[o]=1;
    }else{
        int m=l+r>>1;
        build(ls,l,m);
        build(rs,m+1,r);
        len[o]=len[ls]+len[rs];
        pu(o);
    }
}
void add(int o,int l,int r){
    if(l>R||r<L)return;
    if(l>=L&&r<=R){
        do_1(o,val);
    }else{
        pd(o);
        int m=l+r>>1;
        add(ls,l,m);
        add(rs,m+1,r);
        pu(o);
    }
}
ll get(int o,int l,int r){
    if(l>R||r<L)return 0;
    if(l>=L&&r<=R)return x[o];
    pd(o);
    int m=l+r>>1;
    return (get(ls,l,m)+get(rs,m+1,r))%mod;
}
void add_seg(int l,int r,ll v){
    L=l;R=r;val=v;
    if(L>R)swap(L,R);
    add(1,1,cnt);
}
ll get_seg(int l,int r){
    L=l;R=r;
    if(L>R)swap(L,R);
    return get(1,1,cnt);
}

//operations
void add_path(int x,int y,ll w){
    if(top[x]==top[y])add_seg(id[x],id[y],w);
    else{
        if(deep[top[x]]<deep[top[y]])swap(x,y);
        add_seg(id[top[x]],id[x],w);
        add_path(fa[top[x]],y,w);
    }
}
ll get_path(int x,int y){
    if(top[x]==top[y])return get_seg(id[x],id[y]);
    else{
        if(deep[top[x]]<deep[top[y]])swap(x,y);
        return (get_seg(id[top[x]],id[x])+get_path(fa[top[x]],y))%mod;
    }
}
void add_tree(int x,ll w){
    add_seg(id[x],id[x]+sz[x]-1,w);
}
ll get_tree(int x){
    return get_seg(id[x],id[x]+sz[x]-1);
}

int lca(int x,int y){
    while(top[x]!=top[y]){
        if(deep[top[x]]<deep[top[y]])swap(x,y);
        x=fa[top[x]];
    }
    if(deep[x]<deep[y])return x;
    else return y;
}



int main (){
    ll op,x,y,z;
    scanf("%lld %lld %lld %lld",&n,&m,&rt,&mod);
    for(int i=1;i<=n;++i)scanf("%lld",a+i);
    for(int i=1;i<n;++i){
        scanf("%lld %lld",&x,&y);
        g[x].push_back(y);
        g[y].push_back(x);
    }
    dfs1(rt,0,0);
    dfs2(rt,rt);
    build(1,1,cnt);

    //4种操作
    for(int i=1;i<=m;++i){
        scanf("%lld",&op);
        if(op==1){
            scanf("%lld %lld %lld",&x,&y,&z);
            add_path(x,y,z);
        }else if(op==2){
            scanf("%lld %lld",&x,&y);
            printf("%lld\n",get_path(x,y));
        }else if(op==3){
            scanf("%lld %lld",&x,&y);
            add_tree(x,y);
        }else if(op==4){
            scanf("%lld",&x);
            printf("%lld\n",get_tree(x));
        }
    }
}
```



### Kruskal重构树

```c++
//树上两点的lca的点权即为路径上最大边权
int root[M];
int find(int o){
    return root[o]==o?o:(root[o]=find(root[o]));
}

int val[M*2],cnt;
vector<int> g[M*2];
void kru(){
    cnt=n;
    for(int i=1;i<=n;++i)root[i]=i;
    sort(a+1,a+1+m);//a[M] 边集
    for(int i=1;i<=m;++i)
    {
        int fu=find(a[i].x),fv=find(a[i].y);
        if(fu!=fv)
        {
            val[++cnt]=a[i].w;
            root[cnt]=root[fu]=root[fv]=cnt;
            g[fu].push_back(cnt);g[cnt].push_back(fu);
            g[fv].push_back(cnt);g[cnt].push_back(fv);
        }
    }
}
```



### Dfs序建可持久化线段树

```c++
//https://ac.nowcoder.com/acm/contest/35/D
/*
数据范围: 2e5*1e6
无离散化,值域1e12建树,数组至少开50倍,内存消耗约200M，运行时间牛客600ms. 
*/
#include<bits/stdc++.h>
using namespace std;
const int M=200010;
typedef long long ll;
const ll inf=1e12;
int n,q;
struct nd{
    int x,w;
};
vector<nd> g[M];

ll p,L,R,ans_num;
struct Seg{
    static const int N=M<<6;
    int rt[M],ls[N],rs[N],num[N],cnt;
    ll x[N];

    void pu(int o){
        x[o]=x[ls[o]]+x[rs[o]];
        num[o]=num[ls[o]]+num[rs[o]];
    }
    int clone(int o){
        ++cnt;
        ls[cnt]=ls[o];
        rs[cnt]=rs[o];
        x[cnt]=x[o];
        num[cnt]=num[o];
        return cnt;
    }
    void add(int& o,ll l,ll r){
        o=clone(o);
        if(l==r){
            x[o]+=l;
            num[o]++;
        }else{
            ll m=l+r>>1;
            if(p<=m)add(ls[o],l,m);
            else add(rs[o],m+1,r);
            pu(o);
        }
    }
    ll get(int o,ll l,ll r){
        if(l>R||r<L||o==0)return 0;
        if(l>=L&&r<=R){
            ans_num+=num[o];
            return x[o];
        }else{
            ll m=l+r>>1;
            return get(ls[o],l,m)+get(rs[o],m+1,r);
        }
    }

    void add(int o,ll p0){
        if(p0==0)return;
        p=p0;
        add(rt[o]=rt[o-1],1,inf);
    }
    ll get(int lo,int ro,ll l,ll r){
        L=l;R=r;
        ans_num=0;
        ll x1=get(rt[lo-1],1,inf);
        ans_num*=-1;
        ll x2=get(rt[ro],1,inf);
        return x2-x1;
    }
}seg;

int dfn[M],sz[M],cnt;
ll deep[M];
void dfs1(int o,int fa,ll d){
    dfn[o]=++cnt;
    seg.add(cnt,d);
    deep[o]=d;
    sz[o]=1;
    for(auto v:g[o]){
        if(fa==v.x)continue;
        dfs1(v.x,o,d+v.w);
        sz[o]+=sz[v.x];
    }
}

int main(){
    scanf("%d",&n);
    int x,y;
    for(int i=2;i<=n;++i){
        scanf("%d %d",&x,&y);
        g[i].push_back({x,y});
        g[x].push_back({i,y});
    }
    dfs1(1,0,0);
    scanf("%d",&q);
    for(int i=1;i<=q;++i){
        scanf("%d %d",&x,&y);
        ll t=seg.get(dfn[x],dfn[x]+sz[x]-1,y+deep[x],inf);
        printf("%lld\n",t-ans_num*deep[x]);
    }
}
```





## 图论

### Kosaraju强连通分量

```c++
#include<bits/stdc++.h>
#define M 1005
using namespace std;
int vis[M],sccno[N],scc_cnt;
vector<int> G[M],G2[M];
vector<int> S;

void dfs1(int u){
    if(vis[u])return;
    vis[u]=1;
    for(int i=0;i<G[u].size();++i)dfs1(G[u][i]);
    S.push_back(u);
}

void dfs2(int u){
    if(sccno[u])return;
    sccon[u]=scc_cnt;
    for(int i=0;i<G2[u].size();++i)dfs2(G2[u][i]);
}

void find_scc(int n){
    scc_cnt=0;
    S.clear();
    memset(sccno,0,sizeof(sccno));
    memset(vis,0,sizeof(vis));
    for(int i=0;i<n;++i)dfs1(i);
    for(int i=n-1;i>=0;--i)scc_cnt++,dfs2(S[i]);
}
```



### Tarjan缩点

```c++
#include<bits/stdc++.h>
#define M 1005
using namespace std;
int pre[M],lowlink[M],sccno[M],scc_cnt,dfs_clock,
stack<int> S;
vector<int> G[M];

void dfs(int u){
    pre[u]=lowlink[u]=++dfs_clock;
    S.push(u);
    for(int i=0;i<G[u].size();++i){
        int v=G[u][i];
        if(!pre[v]){
            dfs(v);
            lowlink[u]=min(lowlink[u],lowlink[v]);
        }else if(!sccno[v])//notice..
            lowlink[u]=min(lowlink[u],pre[v]);
    }
    if(lowlink[u]==pre[u]){
        scc_cnt++;
        while(1){
            int x=S.top();S.pop();
            sccno[x]=scc_cnt;
            if(x==u)break;
        }
    }
}

void find_scc(int n){
    dfs_clock=0;
    memset(sccno,0,sizeof(sccno));
    memset(pre,0,sizeof(pre));
    for(int i=1;i<=n;++i)if(!pre[i])dfs(i);
}
```



### 单点k短路

```c++
struct B{
	ll x,t;
	B(ll x,ll t):x(x),t(t){}
	friend bool operator < (B aa,B bb){
		return aa.t+dis[aa.x]>bb.t+dis[bb.x];
		//优先队列，dis为预处理以终点为起点最短路
	}
};

ll bfs(){
	if(dis[s]>=inf)return -1;
	priority_queue <B> q;
	B re(s,0);
	int num=0;
	q.push(re);
	while(!q.empty()){
		re=q.top();
		q.pop();
		if(re.x==e){
			num++;
			if(num==k)return re.t;
		}
		for(int i=pre2[re.x];~i;i=a2[i].next)
			q.push(B(a2[i].y,a2[i].len+re.t));
	}
	return -1;
}
```



### 双连通分量

```
在一个无向图中，若任意两点间至少存在两条“点不重复”的路径，则说这个图是点双连通的。
在一个无向图中，点双连通的极大子图称为点双连通分量。
```

#### 点双连通

```c++
//通用
#include<vector>
#include<string.h>
#include<stdio.h>
using namespace std;
int n,m;

const int maxn=10010;
vector<int> edge[maxn];
vector<vector<int>> connect;
int dfn[maxn],low[maxn],in_seq[maxn];
int stack[maxn],list[maxn];
int cnt,top,pop,len;
void biconnect(int v){
    stack[++top]=v;
    dfn[v]=low[v]=pop++;
    int i,succ;
    for(i=edge[v].size()-1;i>=0;--i){
        succ=edge[v][i];
        if(dfn[succ]==-1){
            biconnect(succ);
            if(low[succ]>=dfn[v]){
                cnt++;
                len=0;
                do{
                    in_seq[stack[top]]=cnt;
                    list[len++]=stack[top];
                    --top;
                }while(stack[top+1]!=succ);
                in_seq[v]=cnt;
                list[len++]=v;
                vector<int> tmp(list,list+len);
                connect.push_back(tmp);
            }
            low[v]=min(low[v],low[succ]);
        }else low[v]=min(low[v],dfn[succ]);
    }
}
void find_bcc(){
    cnt=top=len=pop=0;
    memset(dfn,-1,sizeof dfn);
    connect.clear();
    for(int i=1;i<=n;++i)
        if(dfn[i]==-1)biconnect(i);
}

int main(){
    scanf("%d%d",&n,&m);
    int x,y;
    for(int i=0;i<m;++i){
        scanf("%d%d",&x,&y);
        edge[x].push_back(y);
        edge[y].push_back(x);
    }
    find_bcc();
    for(auto vc:connect){
        for(v:vc)printf("%d ",v);
        puts("");
    }//给定无向图, 输出其所有双连通块.
}

```

```c++
//链式前向星
struct Edge{
    int u,v;
    Edge(int u=0,int v=0):u(u),v(v){}
}e[maxm];
int n,m,stamp,dfn[maxn],low[maxn],iscut[maxn],bccno[maxn];
int scnt,stack[maxm],bcc_cnt;
vector<int> vec[maxn],bcc[maxn];
 
void tarjan(int index,int fa){
    int child=0,tmp;
    dfn[index]=low[index]=++stamp;
    for(int i=0;i<vec[index].size();i++){
        tmp=e[vec[index][i]].v;
        if(!dfn[tmp]){
            stack[++scnt]=vec[index][i],child++;
            tarjan(tmp,index);
            low[index]=min(low[index],low[tmp]);
            if(low[tmp]>=dfn[index]){
                iscut[index]=1;
                bcc[++bcc_cnt].clear();
                while(1){
                    int num=stack[scnt--];
                    if(bccno[e[num].u]!=bcc_cnt){
                        bcc[bcc_cnt].push_back(e[num].u);
                        bccno[e[num].u]=bcc_cnt;
                    }
                    if(bccno[e[num].v]!=bcc_cnt){
                        bcc[bcc_cnt].push_back(e[num].v);
                        bccno[e[num].v]=bcc_cnt;
                    }
                    if(e[num].u==index && e[num].v==tmp)
                        break;
                }
            }
        }
        else if(dfn[tmp]<dfn[index] && tmp!=fa){
            stack[++scnt]=vec[index][i];
            low[index]=min(low[index], dfn[tmp]);
        }
    }
    if(fa<0 && child==1)iscut[index]=0;
}
 
void find_bcc(){
    // 割顶的bccno值无意义 
    memset(dfn,0,sizeof(dfn));
    memset(low,0,sizeof(low));
    memset(iscut,0,sizeof(iscut));
    memset(bccno,0,sizeof(bccno));
    memset(bcc,0,sizeof(bcc));
    stamp=scnt=bcc_cnt=0;
    for(int i=1;i<=n;i++)
        if(!dfn[i])tarjan(i,-1);
}
```

#### 边双连通

```c++
struct Edge{
    int u,v;
    Edge(int u=0,int v=0):u(u),v(v){}
}e[maxm];
int n,m,stamp,dfn[maxn],low[maxn],bccno[maxn],bcc_cnt;
vector<int> vec[maxn],bcc[maxn];
bool g[maxn][maxn],isbridge[maxm];
 
void tarjan(int index,int fa){
    int tmp;
    dfn[index]=low[index]=++stamp;
    for(int i=0;i<vec[index].size();i++){
        tmp=e[vec[index][i]].v;
        if(!dfn[tmp]){
            tarjan(tmp,index);
            low[index]=min(low[index],low[tmp]);
            if(low[tmp]>dfn[index])
                isbridge[vec[index][i]]=isbridge[vec[index][i]^1]=1;
        }
        else if(dfn[tmp]<dfn[index] && tmp!=fa){
            low[index]=min(low[index], dfn[tmp]);
        }
    }
}
void dfs(int index){
    dfn[index]=1;
    bccno[index]=bcc_cnt;
    for(int i=0;i<vec[index].size();i++){
        int tmp=vec[index][i];
        if(isbridge[tmp])continue;
        if(!dfn[e[tmp].v]){
            dfs(e[tmp].v);
        }
    }
}
 
void find_ebcc(){
    bcc_cnt=stamp=0;
    memset(dfn,0,sizeof(dfn));
    memset(low,0,sizeof(low));
    memset(isbridge,0,sizeof(isbridge));
    memset(bccno,0,sizeof(bccno));
    memset(bcc,0,sizeof(bcc));
    for(int i=1;i<=n;i++)
        if(!dfn[i])tarjan(i, -1);
    memset(dfn,0,sizeof(dfn));
    for(int i=1;i<=n;i++){
        if(!dfn[i]){
            bcc_cnt++;
            dfs(i);
        }
    }               
}
```



### 匈牙利算法

```c++
const int M=555;
const int n=100;

vector<int> g[M];
int from[M],tot;
bool use[M];

bool match(int u){
    for(int i=0;i<g[u].size();++i){
        int v=g[u][i];
        if(use[v]==0){
            use[v]=1;
            if(from[v]==-1||match(from[v])){
                from[v]=u;
                return 1;
            }
        }
    }
    return 0;
}

int hungary(){
    tot=0;
    memset(from,255,sizeof from);
    for(int i=1;i<=n;++i){
        memset(use,0,sizeof use);
        if(match(i))++tot;
    }
    return tot;
}
```



### 最短路算法

#### bellman-ford

```c++
bool bellman_ford(int s){
    queue<int> q;
    memset(inq,0,sizeof(inq));
    memset(cnt,0,sizeof(cnt));
    for(int i=0;i<n;++i)d[i]=INF;
    d[s]=0;
    inq[s]=1;
    q.push(s);
    while(!q.empty()){
        int u=q.front();
        q.pop();
        inq[u]=0;
        for(int i=0;i<G[u].size();++i){
            Edge& e=edges[G[u][i]];
            if(d[u]<INF&&d[e.to]>d[u]+e.dist){
                d[e.to]=d[u]+e.dist;
                p[e.to]=G[u][i];
                if(!inq[e.to]){
                    q.push(e.to);
                    inq[e.to]=1;
                    if(++cnt[e.to]>n)return false;
                }
            }
        }
    }
    return true;
}
```



### 网络流

#### Dinic

```c++
//弧优化/自动反向边/封装
const int inf=1<<30;
struct Graph{
    static const int M=100010,N=10010;
    int dis[N],cur[N],head[N],cnt=-1;
    Graph(){
        memset(head,-1,sizeof head);
    }
    struct edge{
        int to,val,next;
    }e[M<<1];
    void _add(int u,int v,int val){
        cnt++;
        e[cnt].next=head[u];
        e[cnt].to=v;
        e[cnt].val=val;
        head[u]=cnt;
    }
    void add(int u,int v,int val){
        _add(u,v,val);
        _add(v,u,0);
    }
    bool bfs(int s,int t){
        queue<int> q;
        memset(dis,-1,sizeof dis);
        dis[s]=0;
        q.push(s);
        while(!q.empty()){
            int u=q.front();
            q.pop();
            for(int i=head[u];i!=-1;i=e[i].next){
                int v=e[i].to;
                if(dis[v]==-1&&e[i].val>0){
                    dis[v]=dis[u]+1;
                    q.push(v);
                }
            }
        }
        return dis[t]!=-1;
    }
    int dfs(int s,int t,int maxflow){
        if(s==t)return maxflow;
        int res=0;
        for(int& i=cur[s];i!=-1;i=e[i].next){
            int v=e[i].to;
            if(dis[v]!=dis[s]+1||e[i].val<=0||res>=maxflow)continue;
            int f=dfs(v,t,min(e[i].val,maxflow-res));
            e[i].val-=f;
            e[i^1].val+=f;
            res+=f;
        }
        return res;
    }
    int Dinic(int s,int t){
        int ans=0;
        while(bfs(s,t)){
            memcpy(cur,head,sizeof head);
            ans+=dfs(s,t,inf);
        }
        return ans;
    }
}mode;
```

```c++
//网络板子
#include<iostream>
#include<cstdio>
#include<cstring>
#include<queue>
#define N 10010
#define M 100010
#define INF 0x7fffffff
using namespace std;
int n,m,ss,tt;
int dis[N],cur[N];
queue<int> q;

struct Edge{int to;int value;int next;}e[M<<1];
int head[N],cnt=-1;
void add(int from,int to,int value){
    cnt++;
    e[cnt].to=to;
    e[cnt].value=value;
    e[cnt].next=head[from];
    head[from]=cnt;
}
bool bfs(int s,int t){//和EK算法的相似,不同的是Dinic中的bfs要求出所有点到源点s的最短路dis[i]
    q=queue<int>();//清空队列
    memset(dis,-1,sizeof(dis));
    dis[s]=0;
    q.push(s);
    while(!q.empty()){
        int x=q.front();
        q.pop();
        for(int i=head[x];i>-1;i=e[i].next){
            int now=e[i].to;
            if(dis[now]==-1&&e[i].value!=0){
                dis[now]=dis[x]+1;
                q.push(now);
            }
        }
    }
    return dis[t]!=-1;
}
int dfs(int x,int t,int maxflow){
    if(x==t)return maxflow;
    int ans=0;
    for(int i=cur[x];i>-1;i=e[i].next){
        int now=e[i].to;
        if(dis[now]!=dis[x]+1||e[i].value==0||ans>=maxflow)continue;
        cur[x]=i;//此路可行，记录此路
        int f=dfs(now,t,min(e[i].value,maxflow-ans));
        e[i].value-=f;
        e[i^1].value+=f;
        ans+=f;
    }
    return ans;
}
int Dinic(int s,int t){
    int ans=0;
    while(bfs(s,t)){
        memcpy(cur,head,sizeof(head));//初始化
        ans+=dfs(s,t,INF);
    }
    return ans;
}

int main(){
    memset(head,-1,sizeof(head));
    scanf("%d%d%d%d",&n,&m,&ss,&tt);
    for(int i=1;i<=m;i++){
        int a,b,c;
        scanf("%d%d%d",&a,&b,&c);
        add(a,b,c);
        add(b,a,0);
    }
    printf("%d\n",Dinic(ss,tt));
    return 0;
}
```





## 字符串

### KMP

```c++
int nxt[M];
void getNext(char b[],int m){
    int i=0,j=-1;
    nxt[0]=-1;
    while(i<m){
        if(j==-1||b[j]==b[i]) nxt[++i]=++j;
        else j=nxt[j];
    }
}
int Kmp(char a[],char b[],int n,int m){
    getNext();
    int i=0,j=0,ans=0;
    while(i<n&&j<m){
        if(j==-1||a[i]==b[j]){
            ++i;++j;
            if(j==m)++ans,j=nxt[j];
        }else j=nxt[j];
    }
    return ans;
}
```





## 其他论题

### 入门经典数位dp

```c++
#include<bits/stdc++.h>
int T,a,b,c,f[35],g[35],h[35];
typedef long long ll;
ll dp[35][2][2][2][2][3][3];
ll Get(int n,int lma,int lmb,int zra,int zrb,int And,int Xor){
    ll &get=dp[n][lma][lmb][zra][zrb][And+1][Xor+1];
    if(get!=-1)return get;
    if(n==0)return get=zra&&zrb&&(And>0||Xor<0);
    get=0;n--;
    for(int i=0;i<=(lma?1:f[n]);++i)for(int j=0;j<=(lmb?1:g[n]);++j)
        get+=Get(n,lma|(i<f[n]),lmb|(j<g[n]),zra|i,zrb|j,And?And:(i&j)-h[n],Xor?Xor:(i^j)-h[n]);
    return get;
}
int main(){
    scanf("%d",&T);
    while(T--){
        memset(dp,-1,sizeof dp);
        scanf("%d%d%d",&a,&b,&c);
        for(int i=0;i<31;++i)f[i]=bool(a&(1<<i)),g[i]=bool(b&(1<<i)),h[i]=bool(c&(1<<i));
        printf("%lld\n",Get(31,0,0,0,0,0,0));
    }
}
```



### LCS-LIS

```c++
//LCS:O(n^2)
int dp[M+1][M+1];
int LCS(int n1,int n2,int A[],int B[]){
    for(int i=1;i<=n1;++i)
    for(int j=1;j<=n2;++j){
        dp[i][j]=dp[i-1][j];
        if(dp[i][j-1]>dp[i][j])dp[i][j]=dp[i][j-1];
        if(A[i]==B[j]&&dp[i-1][j-1]+1>dp[i][j])dp[i][j]=dp[i-1][j-1]+1;
    }
    return dp[n1][n2];
}

//LCIS:O(n^2)
int n1,n2;
int A[M+1],B[M+1];
int dp[M+1][M+1],pre[M+1][M+1],lcis[M+1];
void getLcis(){
    memset(dp,0,sizeof dp);
    memser(pre,0 sizeof pre);
    for(int i=1;i<=n1;++i){
        int k=0;
        for(int j=1;j<=n2;++j){
            if(A[i]!=B[j])dp[i][j]=dp[i-1][j];
            if(A[i]>B[j]&&dp[i][j]>dp[i][k])k=j;
            if(A[i]==B[j]){
                dp[i][j]=dp[i][k]+1;
                pre[i][j]=k;
            }
        }
    }
    int ans=-1;x=n1;y=0;
    for(int i=1;i<=n2;++i)
        if(dp[n1][i]>ans){
            ans=dp[n1][i];
            y=i;
        }
    int cnt=1;
    while(dp[x][y]){
        if(A[x]!=B[y])x--;
        else{
            lcis[ans-cnt]=B[y];
            cnt++;
            y=pre[x][y];
        }
    }
}

//LIS:O(nlongn)
int n,a[M+1],dp[M+1],lis[M+1];
int Lis(){
    memset(dp,88,sizeof dp);
    int ans=0;
    dp[0]=0;
    for(int i=1;i<=n;++i){
        lis[i]=lower_bound(dp+1,dp+n+1,a[i])-dp;
        dp[lis[i]]=a[i];
        ans=max(ans,lis[i]);
    }
    return ans;
}
```



### 模拟退火

```c++
///把n个物品分成m堆，求m个权值和的最小均方差。
#include<bits/stdc++.h>
using namespace std;
int n,m,a[24],s[24];
double T0,delta,Tk,ans,f[24][24];
double sqr(double x){return x*x;}

double work(){//dp跑有序分堆答案
    memset(f,88,sizeof f);
    for(int i=1;i<=n;++i)s[i]=s[i-1]+a[i];
    f[0][0]=0;
    for(int i=1;i<=n;++i)
        for(int j=1;j<=i;j++)
            for(int k=0;k<i;++k)
                f[i][j]=min(f[i][j],f[k][j-1]+sqr(s[i]-s[k]-1.0*s[n]/m));
    return f[n][m];
}
void SA(){//模拟退火随机改变有序顺序
    double t=T0,nowans=ans;
    while(t>Tk){
        int x=0,y=0;
        while(x==y)x=rand()%n+1,y=rand()%n+1;
        swap(a[x],a[y]);
        double now=work();
        double Delta=now-ans;
        if(Delta<0)ans=nowans=now;
        else if(exp(-Delta/t)*RAND_MAX>rand())nowans=now;
        else swap(a[x],a[y]);
        t*=delta;
    }
    for(int i=1;i<=n;++i)
    for(int j=1;j<i;++j){
        swap(a[i],a[j]);
        ans=min(ans,work());
        swap(a[i],a[j]);
    }
}
double Time(){//返回运行时间
    return (double)clock()/CLOCKS_PER_SEC;
}

int main(){
    srand(233);
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    T0=2000;Tk=1e-6;delta=0.993;///--参数--
    ans=1e10;
    while(Time()<0.8)SA();
    printf("%.2f\n",sqrt(ans/m));
}
```

```c++
///上题的[随机化+贪心]解法
///把n个物品分成m堆，求m个权值和的最小均方差。
#include<bits/stdc++.h>
using namespace std;
int n,m,a[24],x[24],sum;
double ans,ave;
double sqr(double x){return x*x;}

void work(){
    memset(x,0,sizeof x);
    for(int i=1;i<=n;++i){
        int p=1;
        for(int j=1;j<=m;++j)if(x[j]<x[p])p=j;
        x[p]+=a[i];
    }
    double now=0;
    for(int i=1;i<=m;++i)now+=sqr(x[i]-ave);
    ans=min(ans,now);
}

int main(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i),sum+=a[i];
    ans=1e10;
    ave=(double)sum/m;
    for(int i=0;i<500000;++i){
        random_shuffle(a+1,a+1+n);
        work();
    }
    printf("%.2f\n",sqrt(ans/m));
}
```









## 数学

### 矩阵快速幂

```c++
#include<bits/stdc++.h>
using namespace std;
const int M=1000010;
typedef long long ll;
int mod,x0,x1,a,b;
char s[M];

//二进制快速幂
template<class T> T Bpow(T a,ll b){
    T res=1;
    while(b){
        if(b&1)res=(res*a)%mod;
        b>>=1;
        a=(a*a)%mod;
    }
    return res;
}
///十进制快速幂
template<class T> T Spow(T a,char s[]){
    T res=1;
    int p=strlen(s)-1;
    while(p>=0){
        ll k=s[p]-'0';
        res=(res*Bpow(a,k))%mod;
        --p;
        a=Bpow(a,10);
    }
    return res;
}
///矩阵
struct matrix{
    int size;
    ll data[8][8];
    void set(int k){size=k;}
    matrix(){
        size=2;
        memset(data,0,sizeof data);
    }
    matrix(int t){
        size=2;
        memset(data,0,sizeof data);
        for(int i=0;i<size;++i)data[i][i]=t;
    }
    matrix operator * (const matrix& t){
        matrix tmp;
        for(int i=0;i<size;++i)
            for(int j=0;j<size;++j)
                for(int k=0;k<size;++k)
                    tmp.data[i][j]+=data[i][k]*t.data[k][j];
        return tmp;
    }
    matrix operator % (ll mod){
        for(int i=0;i<size;++i)for(int j=0;j<size;++j)data[i][j]%=mod;
        return *this;
    }
    void print(){
        for(int i=0;i<size;++i){
        for(int j=0;j<size;++j)
            cout<<data[j][i]<<' ';
            cout<<endl;
        }
    }
};

/// 牛客多校2019-5-B:
int main(){
    scanf("%d%d%d%d",&x0,&x1,&a,&b);
    scanf("%s%d",s,&mod);
    matrix m;
    m.data[0][0]=a;
    m.data[1][0]=b;
    m.data[0][1]=1;
    matrix ans=Spow(m,s);
    //ans.print();
    printf("%d\n",(ans.data[0][1]*x1%mod+ans.data[1][1]*x0%mod)%mod);
}
```



### 组合数

```c++
#include<bits/stdc++.h>
#define M 200010
#define mod 998244353
typedef long long ll;
using namespace std;
//o(n)预处理，o(1)查询
ll f[M],inv[M];
ll qpow(ll a,ll b){
    ll res=1;
    while(b){
        if(b&1)res=res*a%mod;
        a=a*a%mod;b>>=1;
    }
    return res;
}
ll C(int n,int m){
    if(n<0||m<0||m>n)return 0;
    if(m==0||m==n)return 1;
    return f[n]*inv[n-m]%mod*inv[m]%mod;
}
void CInit(){
    f[1]=1;
    for(int i=2;i<M;i++)f[i]=(f[i-1]*i)%mod;
    inv[M-1]=qpow(f[M-1],mod-2);
    for(int i=M-2;i>=1;i--)inv[i]=(inv[i+1]*(i+1))%mod;
}
```

```c++
//递推预处理 o(n^2)预处理，o(1)查询
int C[M][M];
void init(){
    C[0][0]=1;
    for(int i=1;i<M;++i){
        C[i][0]=1;
        for(int j=1;j<=i;++j)
            C[i][j]=(C[i-1][j]+C[i-1][j-1])%mod;
    }
}
```



### 素数筛

```c++
int pri[M],top,vis[M];
//近乎2倍常数的埃筛。。
void prime(){
    for(int i=2;i<M;++i){
        if(vis[i]==0){
            pri[++top]=i;
            for(ll j=1ll*i*i;j<M;j+=i)vis[j]=1;
        }
    }
}
//理论o(n)，跑起来比埃筛慢的线筛。。
void prime_line(){
    for(int i=2;i<M;++i){
        if(vis[i]==0)pri[++top]=i;
        for(int j=1;j<=top&&1ll*i*pri[j]<M;++j){
            vis[i*pri[j]]=1;
            if(i%pri[j]==0)break;
        }
    }
}
```



### 欧拉函数

```c++
//o(sqrt(n))计算
ll euler(ll n){
    ll res=n;
    for(ll i=2;i*i<=n;i++){
        if(n%i==0) res=res/i*(i-1);
        while(n%i==0) n/=i;
    }
    if(n>1)res=res/n*(n-1);
    return res;
}

//o(nlog(n))预处理
int phi[M];
void phi_init(){
    for(int i=1;i<M;++i)phi[i]=i;
    for(int i=2;i<M;++i){
        if(phi[i]==i){
            for(int j=i;j<M;j+=i)
                phi[j]=phi[j]/i*(i-1);
        }
    }
}

//o(n)线筛+欧拉
int vis[M],phi[M],pri[M],top;
void phi_init(){
  phi[1]=1;
  for(int i=2;i<M;++i){
    if(!vis[i]){
      phi[i]=i-1;
      pri[++top]=i;
    }
    for(int j=1;j<=top;++j){
      if(1ll*i*pri[j]>=M)break;
      vis[i*pri[j]] = 1;
      if(i%pri[j]){
        phi[i*pri[j]]=phi[i]*(pri[j]-1);
      }else{
        phi[i*pri[j]]=phi[i]*pri[j];
        break;
      }
    }
  }
}
```



### 扩欧-逆元

```c++
void ex_gcd(ll a,ll b,ll &x,ll &y,ll &d){
    if(!b){d=a,x=1,y=0;}
    else{
        ex_gcd(b,a%b,y,x,d);
        y-=x*(a/b);
    }
}
ll inv(ll t,ll p){
    ll d,x,y;
    ex_gcd(t,p,x,y,d);
    return d==1?(x+p)%p:-1;
}
```



### 线性模方程

```c++
ll solve(ll a,ll b,ll p){
    ll x,y,d;
    ll g=__gcd(a,p);
    if(b%g)return -1;//无解
    b/=g;a/=g;p/=g;
    ex_gcd(a,p,x,y,d);
    return b*(x+p)%p;//返回最小解
}
```



### 中国剩余定理

```c++
//n个方程：x=a[i](mod m[i]) (0<=i<n)
ll china(int n,ll *a,ll *m){
    ll M=1,res=0;
    for(int i=0;i<=n;i++)M*=m[i];
    for(int i=0;i<=n;i++){
        ll w=M/m[i];
        res=(res+w*inv(w,m[i])*a[i])%M;
    }
    return (res+M)%M;
}
//扩展中国剩余定理
ll ex_china(int n,ll *a,ll *m){
    ll res=a[1],r=m[1],x,y,g;
    for(int i=2;i<=n;++i){
        ll t=((a[i]-res)%m[i]+m[i])%m[i];
        ex_gcd(r,m[i],x,y,g);
        if((a[i]-res)%g)return -1;
        x=x*t/g%m[i];//x=mul(x,t/g,m[i]);//防止爆longlong
        res+=r*x;
        r*=m[i]/g;
        res=(res+r)%r;
    }
    return res;
}
```



### 离散对数

```c++
//BSGS
map<int,int>mp;
int bsgs(int a, int b, int p){    //a^x = b (mod P),(a,p)=1，返回x,x>=1; 无解返回-1
    int m=sqrt(p)+1;mp.clear();
    for(register int i=0,res=b;i<m;++i,res=1ll*res*a%p)mp[res]=i;
    for(register int i=1,tmp=qpow(a,m,p),res=tmp;i<=m+1;++i,res=1ll*res*tmp%p)
        if(mp.count(res))return i*m-mp[res];
    return -1;
}
```



### Simpson积分

```c++
double F(double x){
   	///函数F(x)
}
double simpson(double a,double b){
    double c=a+(b-a)/2;
    return (F(a)+4*F(c)+F(b))*(b-a)/6;
}
double asr(double a,double b,double eps,double A){
    double c=a+(b-a)/2;
    double L=simpson(a,c),R=simpson(c,b);
    if(abs(L+R-A)<=15*eps)return L+R+(L+R-A)/15.0;
    return asr(a,c,eps/2,L)+asr(c,b,eps/2,R);
}
double asr(double a,double b,double eps){
    return asr(a,b,eps,simpson(a,b));
}
```





## 计算几何

### 三角形

```c++
//已知三边平方，计算夹角
double calR(double a2,double b2,double c2){
	return acos((a2+b2-c2)/(2*sqrt(a2*b2)));
}

//已知三边，计算底边高
double calH(double a,double b,double c){
    double r=calR(a*a,b*b,c*c);
    return a*b*sin(r)/c;
}
```



### 点线

```c++
//二维坐标旋转
void calC(double& x,double& y,double r){
    double tx=t,ty=y;
    x=cos(r)*tx-sin(r)*ty;
    y=sin(r)*tx+cos(r)*ty;
}
```

