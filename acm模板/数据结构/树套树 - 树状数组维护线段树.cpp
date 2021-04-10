#include<bits/stdc++.h>
using namespace std;
const int M=600010;
int n,m,a[M],rt[M],rk[M],top;
//树状数组维护动态开点线段树 - 动态区间第k小
int x[M<<5],ls[M<<5],rs[M<<5],cnt;
int lowbit(int x){return x&-x;}
int p,L,R,val;
void pu(int& o){x[o]=x[ls[o]]+x[rs[o]];}
void pd(int& o){if(!o)o=++cnt;}
void _add(int& o,int l,int r){
    pd(o);
    if(l==r)x[o]+=val;
    else{
        int m=l+r>>1;
        if(p<=m)_add(ls[o],l,m);
        else _add(rs[o],m+1,r);
        pu(o);
    }
}
void add(int o,int p0,int val0){//修改权值线段树
    p=p0;val=val0;
    for(;o<=top;o+=lowbit(o))_add(rt[o],1,top);
}
int tmp[50][2];
int _getKth(int ln,int rn,int l,int r,int k){
    if(l==r)return rk[l];
    int lx=0,rx=0;
    for(int i=1;i<=ln;++i)lx+=x[ls[tmp[i][0]]];
    for(int i=1;i<=rn;++i)rx+=x[ls[tmp[i][1]]];
    int m=l+r>>1;
    if((rx-lx)>=k){
        for(int i=1;i<=ln;++i)tmp[i][0]=ls[tmp[i][0]];
        for(int i=1;i<=rn;++i)tmp[i][1]=ls[tmp[i][1]];
        return _getKth(ln,rn,l,m,k);
    }else{
        for(int i=1;i<=ln;++i)tmp[i][0]=rs[tmp[i][0]];
        for(int i=1;i<=rn;++i)tmp[i][1]=rs[tmp[i][1]];
        return _getKth(ln,rn,m+1,r,k-(rx-lx));
    }
}
int getKth(int lo,int ro,int k){//lo-ro区间第k小
    int ln=0,rn=0;
    for(;lo>0;lo-=lowbit(lo))tmp[++ln][0]=rt[lo];
    for(;ro>0;ro-=lowbit(ro))tmp[++rn][1]=rt[ro];
    return _getKth(ln,rn,1,top,k);
}

struct nd{
    int op,l,r,k;
}b[M];

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        rk[++top]=a[i];
    }
    char s[5];
    for(int i=1;i<=m;++i){
        scanf("%s",s);
        if(s[0]=='Q'){
            b[i].op=0;
            scanf("%d %d %d",&b[i].l,&b[i].r,&b[i].k);
        }else{
            b[i].op=1;
            scanf("%d %d",&b[i].l,&b[i].r);
            rk[++top]=b[i].r;
        }
    }
    sort(rk+1,rk+1+top);
    top=unique(rk+1,rk+1+top)-rk-1;
    for(int i=1;i<=n;++i){
        a[i]=lower_bound(rk+1,rk+1+top,a[i])-rk;
        add(i,a[i],1);
    }

    for(int i=1;i<=m;++i){
        if(b[i].op==0){
            printf("%d\n",getKth(b[i].l-1,b[i].r,b[i].k));
        }else{
            int t1=a[b[i].l];
            int t2=lower_bound(rk+1,rk+1+top,b[i].r)-rk;
            a[b[i].l]=t2;
            add(b[i].l,t1,-1);
            add(b[i].l,t2,1);
        }
    }
}
