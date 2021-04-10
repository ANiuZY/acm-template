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
