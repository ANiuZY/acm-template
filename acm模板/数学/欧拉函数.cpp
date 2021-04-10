//o(sqrt(n))计算
ll euler(ll n){
    ll res=n;
    for(ll i=2;i*i<=n;i++){
        if(n%i==0) res=res/i*(i-1);
        while(n%i==0) n/=i;
    }
    if(n>1)res=res/n*(n-1);
    return res;
}

//o(nlog(n))预处理
int phi[M];
void phi_init(){
    for(int i=1;i<M;++i)phi[i]=i;
    for(int i=2;i<M;++i){
        if(phi[i]==i){
            for(int j=i;j<M;j+=i)
                phi[j]=phi[j]/i*(i-1);
        }
    }
}

//o(n)线筛+欧拉
int vis[M],phi[M],pri[M],top;
void phi_init(){
  phi[1]=1;
  for(int i=2;i<M;++i){
    if(!vis[i]){
      phi[i]=i-1;
      pri[++top]=i;
    }
    for(int j=1;j<=top;++j){
      if(1ll*i*pri[j]>=M)break;
      vis[i*pri[j]] = 1;
      if(i%pri[j]){
        phi[i*pri[j]]=phi[i]*(pri[j]-1);
      }else{
        phi[i*pri[j]]=phi[i]*pri[j];
        break;
      }
    }
  }
}
