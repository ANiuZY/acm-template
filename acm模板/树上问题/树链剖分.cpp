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
