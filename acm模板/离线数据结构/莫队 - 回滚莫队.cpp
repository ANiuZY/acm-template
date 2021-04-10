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
