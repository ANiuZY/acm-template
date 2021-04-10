//o(1)子树剖离技巧

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
