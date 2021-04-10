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
