#include<bits/stdc++.h>
using namespace std;
const int M=150010;
int n,m;

int rt[M<<1];
struct seg{
    int x[M<<4],ls[M<<4],rs[M<<4],cnt;
    int deep[M<<4];
    int dp;
    void build(int& o,int l,int r){
        o=++cnt;
        if(l==r){
            x[o]=l;
            deep[o]=1;
        }else{
            int m=l+r>>1;
            build(ls[o],l,m);
            build(rs[o],m+1,r);
        }
    }
    void init(){
        cnt=0;
        build(rt[cnt++],1,n);
    }
    int clone(int o){
        ++cnt;
        x[cnt]=x[o];
        ls[cnt]=ls[o];
        rs[cnt]=rs[o];
        return cnt;
    }
    void add(int& o,int l,int r,int& p,int& val){
        o=clone(o);
        if(l==r){
            x[o]=val;
        }else{
            int m=l+r>>1;
            if(p<=m)add(ls[o],l,m,p,val);
            else add(rs[o],m+1,r,p,val);
        }
    }
    void addDeep(int o,int l,int r,int p,int val){
        if(l==r){
            deep[o]=val;
        }else{
            int m=l+r>>1;
            if(p<=m)addDeep(ls[o],l,m,p,val);
            else addDeep(rs[o],m+1,r,p,val);
        }
    }
    int get(int o,int l,int r,int& p){
        if(l==r){
            dp=deep[o];
            return x[o];
        }
        int m=l+r>>1;
        if(p<=m)return get(ls[o],l,m,p);
        else return get(rs[o],m+1,r,p);
    }
    int Root(int ver,int p){//查询
        int rp=get(rt[ver],1,n,p);
        if(p==rp)return p;
        else return Root(ver,rp);
    }
    void OrderMerge(int ver,int p1,int p2){//按秩合并
        int dp1,dp2,r1,r2;
        r1=Root(ver,p1);dp1=dp;
        r2=Root(ver,p2);dp2=dp;
        if(r1==r2)return;
        if(dp1>dp2)add(rt[ver],1,n,r2,r1);
        else{
            add(rt[ver],1,n,r1,r2);
            if(dp1==dp2)addDeep(rt[ver],1,n,r2,dp2+1);
        }
    }
}S;

int main(){
    scanf("%d %d",&n,&m);
    S.init();
    int op,x1,x2;
    for(int i=1;i<=m;++i){
        scanf("%d",&op);
        if(op==1){
            scanf("%d %d",&x1,&x2);//合并
            rt[i]=rt[i-1];
            S.OrderMerge(i,x1,x2);
        }else if(op==2){
            scanf("%d",&x1);//跳转版本
            rt[i]=rt[x1];
        }else {
            scanf("%d %d",&x1,&x2);//查询
            rt[i]=rt[i-1];
            int r1=S.Root(i,x1);
            int r2=S.Root(i,x2);
            printf("%d\n",r1==r2);
        }
    }
}
