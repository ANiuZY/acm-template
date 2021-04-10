#include<bits/stdc++.h>
const int M=200010;
using namespace std;
typedef long long ll;
ll n,m,a[M];

int rt[M],cnt;
ll x[M<<5];
int ls[M<<5],rs[M<<5];
void pd(int& o){
    if(!o)o=++cnt;
}
void pu(int& o){
    x[o]=x[ls[o]]+x[rs[o]];
}
void build(int& o,int l,int r){
    o=++cnt;
    if(l==r){
        x[o]=a[l];
    }else{
        int m=l+r>>1;
        build(ls[o],l,m);
        build(rs[o],m+1,r);
        pu(o);
    }
}
int L,R,p;
ll val,k;
void move(int& lo,int& ro,int l,int r){//分裂: lo->ro
    if(l>R||r<L)return;
    if(!lo)return;
    if(l>=L&&r<=R){
        swap(ro,lo);
    }else{
        pd(ro);
        int m=l+r>>1;
        move(ls[lo],ls[ro],l,m);
        move(rs[lo],rs[ro],m+1,r);
        pu(ro);
        pu(lo);
    }
}
void merge(int& lo,int& ro){//合并: lo->ro
    if(!lo)return;
    if(!ro)swap(lo,ro);
    else{
        x[ro]+=x[lo];
        x[lo]=0;
        merge(ls[lo],ls[ro]);
        merge(rs[lo],rs[ro]);
    }
}
void add(int& o,int l,int r){
    pd(o);
    if(l==r){
        x[o]+=val;
    }else{
        int m=l+r>>1;
        if(p<=m)add(ls[o],l,m);
        else add(rs[o],m+1,r);
        pu(o);
    }
}
ll get(int& o,int l,int r){
    if(l>R||r<L)return 0;
    pd(o);
    if(l>=L&&r<=R){
        return x[o];
    }else{
        int m=l+r>>1;
        return get(ls[o],l,m)+get(rs[o],m+1,r);
    }
}
int getKth(int& o,int l,int r,ll k){
    pd(o);
    if(l==r){
        return l;
    }else{
        if(x[o]<k)return -1;
        int m=l+r>>1;
        if(k<=x[ls[o]])return getKth(ls[o],l,m,k);
        else return getKth(rs[o],m+1,r,k-x[ls[o]]);
    }
}

int main(){
    scanf("%lld %lld",&n,&m);
    for(int i=1;i<=n;++i)scanf("%lld",a+i);
    int op,x,y;
    int top=1;
    build(rt[top],1,n);
    for(int i=1;i<=m;++i){
        scanf("%d",&op);
        if(op==0){
            scanf("%d %d %d",&x,&L,&R);
            move(rt[x],rt[++top],1,n);//分裂
        }else if(op==1){
            scanf("%d %d",&x,&y);
            merge(rt[y],rt[x]);//合并
        }else if(op==2){
            scanf("%d %lld %lld",&x,&val,&p);
            add(rt[x],1,n);//单点修改
        }else if(op==3){
            scanf("%d %d %d",&x,&L,&R);
            printf("%lld\n",get(rt[x],1,n));//区间查询
        }else {
            scanf("%d %lld",&x,&k);
            printf("%d\n",getKth(rt[x],1,n,k));//第k小
        }
    }
}
