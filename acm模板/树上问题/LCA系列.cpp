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
