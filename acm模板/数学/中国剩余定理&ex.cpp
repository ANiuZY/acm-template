//n个方程：x=a[i](mod m[i]) (0<=i<n)
ll china(int n,ll *a,ll *m){
    ll M=1,res=0;
    for(int i=0;i<=n;i++)M*=m[i];
    for(int i=0;i<=n;i++){
        ll w=M/m[i];
        res=(res+w*inv(w,m[i])*a[i])%M;
    }
    return (res+M)%M;
}

//扩展中国剩余定理
ll ex_china(int n,ll *a,ll *m){
    ll res=a[1],r=m[1],x,y,g;
    for(int i=2;i<=n;++i){
        ll t=((a[i]-res)%m[i]+m[i])%m[i];
        ex_gcd(r,m[i],x,y,g);
        if((a[i]-res)%g)return -1;
        x=x*t/g%m[i];//x=mul(x,t/g,m[i]);//防止爆longlong
        res+=r*x;
        r*=m[i]/g;
        res=(res+r)%r;
    }
    return res;
}
