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
