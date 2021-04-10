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
