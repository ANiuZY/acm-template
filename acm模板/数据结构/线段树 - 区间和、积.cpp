#include<bits/stdc++.h>
using namespace std;
const int M=100010;
typedef long long ll;
ll n,m,p;
ll a[M];

struct nd{
    ll sum,lzs,lzm;
    ll len;
}tr[M<<2];
ll L,R,val;

void pu(ll o){
    tr[o].sum=(tr[o<<1].sum+tr[o<<1|1].sum)%p;
}

void pd(ll o){
    if(tr[o].lzm!=1){
        ll v=tr[o].lzm;
        tr[o<<1].sum=(tr[o<<1].sum*v)%p;
        tr[o<<1|1].sum=(tr[o<<1|1].sum*v)%p;
        tr[o<<1].lzm=(tr[o<<1].lzm*v)%p;
        tr[o<<1].lzs=(tr[o<<1].lzs*v)%p;
        tr[o<<1|1].lzm=(tr[o<<1|1].lzm*v)%p;
        tr[o<<1|1].lzs=(tr[o<<1|1].lzs*v)%p;
        tr[o].lzm=1;
    }
    if(tr[o].lzs){
        ll v=tr[o].lzs;
        tr[o<<1].sum=(tr[o<<1].sum+v*tr[o<<1].len)%p;
        tr[o<<1|1].sum=(tr[o<<1|1].sum+v*tr[o<<1|1].len)%p;
        tr[o<<1].lzs=(tr[o<<1].lzs+v)%p;
        tr[o<<1|1].lzs=(tr[o<<1|1].lzs+v)%p;
        tr[o].lzs=0;
    }
}

void build(ll o,ll l,ll r){
    tr[o].lzm=1;
    if(l==r){
        tr[o].sum=a[l];
        tr[o].len=1;
        return;
    }
    ll m=l+r>>1;
    build(o<<1,l,m);
    build(o<<1|1,m+1,r);
    tr[o].len=tr[o<<1].len+tr[o<<1|1].len;
    pu(o);
}

void add(ll o,ll l,ll r){
    if(r<L||l>R)return;
    if(r<=R&&l>=L){
        tr[o].sum=(tr[o].sum+val*(r-l+1))%p;
        tr[o].lzs=(tr[o].lzs+val)%p;
    }else{
        pd(o);
        ll m=l+r>>1;
        add(o<<1,l,m);
        add(o<<1|1,m+1,r);
        pu(o);
    }
}

void mul(ll o,ll l,ll r){
    if(r<L||l>R)return;
    if(r<=R&&l>=L){
        tr[o].sum=tr[o].sum*val%p;
        tr[o].lzm=tr[o].lzm*val%p;
        tr[o].lzs=tr[o].lzs*val%p;
    }else{
        pd(o);
        ll m=l+r>>1;
        mul(o<<1,l,m);
        mul(o<<1|1,m+1,r);
        pu(o);
    }
}

ll get(ll o,ll l,ll r){
    if(l>R||r<L)return 0;
    if(l>=L&&r<=R)return tr[o].sum;
    pd(o);
    ll m=l+r>>1;
    return (get(o<<1,l,m)+get(o<<1|1,m+1,r))%p;
}

int main(){
    scanf("%lld %lld %lld",&n,&m,&p);
    for(int i=1;i<=n;++i)scanf("%lld",a+i);
    build(1,1,n);
    ll x,y,z;
    for(int i=0;i<m;++i){
        scanf("%lld",&x);
        if(x==1){
            scanf("%lld %lld %lld",&x,&y,&z);
            L=x;R=y;val=z;
            mul(1,1,n);
        }else if(x==2){
            scanf("%lld %lld %lld",&x,&y,&z);
            L=x;R=y;val=z;
            add(1,1,n);
        }else{
            scanf("%lld %lld",&x,&y);
            L=x;R=y;
            printf("%lld\n",get(1,1,n));
        }
    }
}
