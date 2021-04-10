//静态区间众数
//值域1e9,离线求任意区间最多的数的个数
#include<bits/stdc++.h>
using namespace std;
const int M=500010;
const int maxn=1e9+7;
int n,m,a[M],r[M];
struct nd{
    int l,r,L,id;
    bool operator < (const nd& t)const {
        if(L!=t.L)return L<t.L;
        return (r<t.r);
    }
}b[M];

int now,cnt[M],ans[M],num[M];
void add(int t){
    num[cnt[t]]--;
    cnt[t]++;
    num[cnt[t]]++;
    if(cnt[t]>now)now=cnt[t];
}
void del(int t){
    num[cnt[t]]--;
    cnt[t]--;
    num[cnt[t]]++;
    if(cnt[t]+1==now)
        if(num[now]==0)now--;
}

int main(){
    scanf("%d %d" ,&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        r[i]=a[i];
    }
    sort(r+1,r+1+n);
    int nr=unique(r+1,r+1+n)-r-1;
    for(int i=1;i<=n;++i)a[i]=lower_bound(r+1,r+1+nr,a[i])-r;

    int s=sqrt(n);
    for(int i=1;i<=m;++i){
        scanf("%d %d",&b[i].l,&b[i].r);
        b[i].id=i;
        b[i].L=b[i].l/s;
    }
    sort(b+1,b+1+m);
    int l=1,r=0;
    for(int i=1;i<=m;++i){
        int ll=b[i].l,rr=b[i].r;
        while(r<rr)add(a[++r]);
        while(l>ll)add(a[--l]);
        while(r>rr)del(a[r--]);
        while(l<ll)del(a[l++]);
        ans[b[i].id]=now;
    }
    for(int i=1;i<=m;i++)printf("%d\n",ans[i]);
}
