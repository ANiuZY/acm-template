//区间异或v 区间查询任意个数最大异或值
//用差分数组建立线段树
#include<bits/stdc++.h>
#define ls o<<1
#define rs o<<1|1
using namespace std;
const int M=50010;
const int N=31;
int n,m,a[M],b[M];

struct base{
    int a[N];
    base(){memset(a,0,sizeof a);}
    void add(int x){
        for(int i=N-1;i>=0&&x;--i)
            if(x>>i){
                if(!a[i]){
                    a[i]=x;
                    return;
                }
                x^=a[i];
            }
    }
    int get(int x){
        for(int i=N-1;i>=0;--i)x=max(x,x^a[i]);
        return x;
    }
    base merge(base t){
        for(int i=0;i<N;++i)t.add(a[i]);
        return t;
    }
}bs[M<<2];
int sum[M<<2];
int p,val,L,R;
void pu(int o){
    sum[o]=sum[ls]^sum[rs];
    bs[o]=bs[ls].merge(bs[rs]);
}
void build(int o,int l,int r){
    if(l==r){
        sum[o]=b[l];
        bs[o].add(sum[o]);
    }else{
        int m=l+r>>1;
        build(ls,l,m);
        build(rs,m+1,r);
        pu(o);
    }
}
int getSum(int o,int l,int r,int L,int R){
    if(l>R||r<L)return 0;
    if(l>=L&&r<=R)return sum[o];
    int m=l+r>>1;
    return getSum(ls,l,m,L,R)^getSum(rs,m+1,r,L,R);
}
void _add(int o,int l,int r){
    if(l==r){
        sum[o]^=val;
        bs[o]=base();
        bs[o].add(sum[o]);
    }else{
        int m=l+r>>1;
        if(p<=m)_add(ls,l,m);
        else _add(rs,m+1,r);
        pu(o);
    }
}
void add(int l,int r,int v){ //[l,r]区间异或v
    val=v;p=l;
    _add(1,1,n);
    p=r+1;
    if(p<=n)_add(1,1,n);
}
base _get(int o,int l,int r){
    if(l>R||r<L)return base();
    if(l>=L&&r<=R)return bs[o];
    int m=l+r>>1;
    return _get(ls,l,m).merge(_get(rs,m+1,r));
}
int get(int l,int r,int v){ //[l,r]内任意个数与v异或的最大值
    L=l+1;R=r;
    base t;
    if(L<=R)t=_get(1,1,n);
    t.add(getSum(1,1,n,1,l));
    return t.get(v);
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        b[i]=a[i-1]^a[i];//用差分数组建立线段树。
    }
    build(1,1,n);
    int op,l,r,v;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d %d",&op,&l,&r,&v);
        if(op==1) add(l,r,v);
        else printf("%d\n",get(l,r,v));
    }
}
