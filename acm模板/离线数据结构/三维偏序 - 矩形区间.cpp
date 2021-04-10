// 矩形区间三维偏序
#include<bits/stdc++.h>
using namespace std;
const int M=400010;
int n,k,top;
int rk[M];
long long ans;
struct nd{
    int x,r,q;
    int vl,vr,ql,qr;
    bool operator < (const nd& t) const {
        return r>t.r;
    }
}a[M];

bool cmp1(nd& t1,nd& t2){
    return t1.q<t2.q;
}

struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int v){
        for(;x<=top;x+=lowbit(x))
            a[x]+=v;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))
            res+=a[x];
        return res;
    }
}bt;

void cdq(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    cdq(l,m);
    cdq(m+1,r);
    sort(a+l,a+m+1,cmp1);
    sort(a+m+1,a+r+1,cmp1);

    //q的匹配范围为长度为2k的定长区间，所以可以直接分治q;
    //或者将x的不定长匹配区间拆开（左、中、右），分治x。
    //这里分治q
    int pl=l,pr=l,x=m+1;
    while(x<=r){
        while(pr<=m&&a[pr].q<=a[x].q+k)bt.add(a[pr].x, 1),pr++;
        while(pl<=m&&a[pl].q< a[x].q-k)bt.add(a[pl].x,-1),pl++;
        ans+=bt.get(a[x].vr)-bt.get(a[x].vl-1);
        x++;
    }
    while(pl<pr)bt.add(a[pl].x,-1),pl++;
}

int main(){
    scanf("%d %d",&n,&k);
    for(int i=1;i<=n;++i){
        scanf("%d %d %d",&a[i].x,&a[i].r,&a[i].q);
        rk[++top]=a[i].x;
        rk[++top]=a[i].x+a[i].r;
        rk[++top]=a[i].x-a[i].r;
    }
    sort(rk+1,rk+1+top);
    top=unique(rk+1,rk+1+top)-rk-1;
    for(int i=1;i<=n;++i){
        a[i].vl=lower_bound(rk+1,rk+top+1,a[i].x-a[i].r)-rk;
        a[i].vr=lower_bound(rk+1,rk+top+1,a[i].x+a[i].r)-rk;
        a[i].x=lower_bound(rk+1,rk+top+1,a[i].x)-rk;
    }
    sort(a+1,a+1+n);
    cdq(1,n);
    printf("%lld",ans);
}
