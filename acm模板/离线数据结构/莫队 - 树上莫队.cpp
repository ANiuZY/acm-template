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
