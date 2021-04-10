#include<bits/stdc++.h>
using namespace std;
const int N=200010;
const int M=1000010;
int n,m,a[N],ts[N],b[N],ans[N];
int top,cnt;

struct nd{
    int l,r,L,R,t,id;
    bool operator <(const nd& v)const{
        if(L!=v.L)return L<v.L;
        if(R!=v.R)return (R<v.R);
        return (t<v.t);
    }
}q[N];

int now,num[M];
void add(int x){
    if(num[x]==0)now++;
    num[x]++;
}
void del(int x){
    num[x]--;
    if(num[x]==0)now--;
}
void update(int x,int l,int r){
    int i=ts[x];
    if(i>=l&&i<=r){
        add(b[x]);
        del(a[i]);
    }
    swap(a[i],b[x]);
}

int main(){
    char c[10];
    int x,y;
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    for(int i=1;i<=m;++i){
        scanf("%s %d %d",c,&x,&y);
        if(c[0]=='R'){
            ts[++top]=x;
            b[top]=y;
        }else{
            ++cnt;
            q[cnt]=(nd){x,y,0,0,top,cnt};
        }
    }

    int s=pow(n,0.666);
    for(int i=1;i<=cnt;++i){
        q[i].L=q[i].l/s;
        q[i].R=q[i].r/s;
    }
    sort(q+1,q+1+cnt);
    int l=1,r=0,t=0;
    for(int i=1;i<=cnt;++i){
        int x=q[i].l,y=q[i].r,z=q[i].t;
        while(r<y)add(a[++r]);
        while(l>x)add(a[--l]);
        while(r>y)del(a[r--]);
        while(l<x)del(a[l++]);
        while(t<z)update(++t,l,r);
        while(t>z)update(t--,l,r);
        ans[q[i].id]=now;
    }
    for(int i=1;i<=cnt;++i)printf("%d\n",ans[i]);
}
