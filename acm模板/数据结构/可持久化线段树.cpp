#include<bits/stdc++.h>
using namespace std;
const int M=1000010;
int n,m,cnt,a[M],rt[M];

struct nd{
    int l,r,x;
}tr[M<<4];

void build(int& o,int l,int r){
    o=++cnt;
    if(l==r)tr[o].x=a[l];
    else {
        int m=l+r>>1;
        build(tr[o].l,l,m);
        build(tr[o].r,m+1,r);
    }
}
int clone(int o){
    tr[++cnt]=tr[o];
    return cnt;
}
int p,val,ver;
void add(int& o,int l,int r){
    o=clone(o);
    if(l==r){
        tr[o].x=val;
    }else{
        int m=l+r>>1;
        if(p<=m)add(tr[o].l,l,m);
        else add(tr[o].r,m+1,r);
    }
}
int get(int o,int l,int r){
    if(l==r){
        return tr[o].x;
    }else{
        int m=l+r>>1;
        if(p<=m)return get(tr[o].l,l,m);
        else return get(tr[o].r,m+1,r);
    }
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    build(rt[0],1,n);
    int op;
    for(int i=1;i<=m;++i){
        scanf("%d %d",&ver,&op);
        rt[i]=rt[ver];
        if(op==1){
            scanf("%d %d",&p,&val);//修改
            add(rt[i],1,n);
        }else{
            scanf("%d",&p);
            printf("%d\n",get(rt[i],1,n));// 查询
        }
    }
}
