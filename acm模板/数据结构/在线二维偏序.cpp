//在线二维偏序(p4602)
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int M=100010;
const int N=100000;
int n,m;
struct nd{
    int d,p,l;
    bool operator < (const nd& t)const {
        return d<t.d;
    }
}a[M];
struct nd2{
    int l,r;
    ll x,p;
}tr[M*50];

int rt[M],cnt,val,p;
int clone(int o){
    tr[++cnt]=tr[o];
    return cnt;
}
void pu(int o){
    tr[o].x=tr[tr[o].l].x+tr[tr[o].r].x;
    tr[o].p=tr[tr[o].l].p+tr[tr[o].r].p;
}
void add(int &o,int l,int r){
    o=clone(o);
    if(l==r){
        tr[o].x+=a[p].l;
        tr[o].p=tr[o].x*l;
    }else{
        int m=l+r>>1;
        if(a[p].p<=m)add(tr[o].l,l,m);
        else add(tr[o].r,m+1,r);
        pu(o);
    }
}
ll get(int o,int l,int r,ll p){
    if(l==r)return min(p/l,tr[o].x);
    int m=l+r>>1;
    ll t=tr[tr[o].l].p;
    if(t<=p)return tr[tr[o].l].x+get(tr[o].r,m+1,r,p-t);
    else return get(tr[o].l,l,m,p);
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)
        scanf("%d %d %d",&a[i].d,&a[i].p,&a[i].l);
    sort(a+1,a+n+1);
    for(p=n;p>0;--p){
        rt[p]=rt[p+1];
        add(rt[p],1,N);
    }

    ll x,y;
    for(int i=1;i<=m;++i){
        scanf("%lld %lld",&x,&y);
        int l=0,r=n;
        while(l<r){
            int mid=(l+r+1)>>1;
            if(get(rt[mid],1,N,x)>=y)l=mid;
            else r=mid-1;
        }
        if(l)printf("%d\n",a[l].d);
        else printf("-1\n");
    }
}
