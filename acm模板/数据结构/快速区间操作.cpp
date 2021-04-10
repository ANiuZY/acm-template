//log修改 log查询。
struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int val){
        for(;x<=n;x+=lowbit(x))a[x]+=val;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))res+=a[x];
        return res;
    }
}tr;

//sqrt(n)修改，o(1)查询。
struct BT{
    int k[M],as[M],aq[M],s,ns;
    void init(){// 使用前先init()!!!
        s=sqrt(n);ns=n/s+1;
        for(int i=1;i<=n;++i)k[i]=i/s+1;
    }
    void add(int x,int val){
        for(int i=k[x]+1;i<=ns;++i)as[i]+=val;
        for(int x0=x;x<=n&&k[x]==k[x0];++x)aq[x]+=val;
    }
    int get(int x){
        return as[k[x]]+aq[x];
    }
}bt;
