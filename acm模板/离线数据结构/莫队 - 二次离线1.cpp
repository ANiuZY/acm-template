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
