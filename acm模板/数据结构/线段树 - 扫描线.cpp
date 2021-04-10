#include<bits/stdc++.h>
#define ls o<<1
#define rs o<<1|1
using namespace std;
typedef long long ll;
const int M=400010;
int n,rk[M],top,cnt;
struct nd{
    int x,y1,y2,f;
    nd(){}
    nd(int x1,int x2,int x3,int x4){x=x1;y1=x2;y2=x3;f=x4;}
    bool operator < (const nd& t)const {return x<t.x;}
}a[M];
//扫描线线段树
int x[M<<2],num[M<<2],v[M<<2];
void pu(int o){
    if(num[o])x[o]=v[o];
    else x[o]=x[ls]+x[rs];
}
void build(int o,int l,int r){
    if(l==r){
        v[o]=rk[l+1]-rk[l];
    }else{
        int m=l+r>>1;
        build(ls,l,m);
        build(rs,m+1,r);
        v[o]=v[ls]+v[rs];
    }
}
int L,R,val;
void add(int o,int l,int r){
    if(l>R||r<L)return;
    if(l>=L&&r<=R){
        num[o]+=val;
    }else{
        int m=l+r>>1;
        add(ls,l,m);
        add(rs,m+1,r);
    }
    pu(o);// !
}

int main(){
    scanf("%d",&n);
    int x1,x2,y1,y2;
    for(int i=1;i<=n;++i){
        scanf("%d %d %d %d",&x1,&y1,&x2,&y2);
        rk[++top]=y1;
        rk[++top]=y2;
        a[++cnt]=nd(x1,y1,y2,1);
        a[++cnt]=nd(x2,y1,y2,-1);
    }
    sort(a+1,a+1+cnt);
    sort(rk+1,rk+1+top);
    top=unique(rk+1,rk+1+top)-rk-1;
    build(1,1,top-1);
    ll ans=0;
    int p=1;
    for(int i=1;i<=cnt;++i){
        a[i].y1=lower_bound(rk+1,rk+1+top,a[i].y1)-rk;
        a[i].y2=lower_bound(rk+1,rk+1+top,a[i].y2)-rk;
    }
    for(int i=1;i<=cnt;++i){
        ans+=1ll*(a[i].x-a[i-1].x)*x[1];
        L=a[i].y1;R=a[i].y2-1;val=a[i].f;
        add(1,1,top-1);
    }
    printf("%lld\n",ans);
}
