int pri[M],top,vis[M];
//近乎2倍常数的埃筛。。
void prime(){
    for(int i=2;i<M;++i){
        if(vis[i]==0){
            pri[++top]=i;
            for(ll j=1ll*i*i;j<M;j+=i)vis[j]=1;
        }
    }
}
//理论o(n)，跑起来比埃筛慢的线筛。。
void prime_line(){
    for(int i=2;i<M;++i){
        if(vis[i]==0)pri[++top]=i;
        for(int j=1;j<=top&&1ll*i*pri[j]<M;++j){
            vis[i*pri[j]]=1;
            if(i%pri[j]==0)break;
        }
    }
}
