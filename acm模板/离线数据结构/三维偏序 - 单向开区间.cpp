//cdq分治模板题 单向开区间三维偏序
#include<bits/stdc++.h>
using namespace std;
const int M=200010;
int n,top,m,ans[M];

struct nd{
    int x,y,z,cnt,ans;
}a[M],b[M];
bool cmpy(const nd& t1,const nd& t2){
    if(t1.y==t2.y)return t1.z<t2.z;
    return t1.y<t2.y;
}
bool cmpx(const nd& t1,const nd& t2){
    if(t1.x==t2.x)return cmpy(t1,t2);
    return t1.x<t2.x;
}

struct BIT{
    int a[M];
    int lowbit(int x){return x&-x;}
    void add(int x,int v){
        for(;x<=m;x+=lowbit(x))a[x]+=v;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))res+=a[x];
        return res;
    }
}bt;

void cdq(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    cdq(l,m);
    cdq(m+1,r);
    sort(a+l,a+m+1,cmpy);
    sort(a+m+1,a+r+1,cmpy);

    int x=l,y=m+1;
    while(y<=r){
        if(x<=m&&a[x].y<=a[y].y){
            bt.add(a[x].z,a[x].cnt);
            x++;
        }else{
            a[y].ans+=bt.get(a[y].z);
            y++;
        }
    }
    for(;l<x;++l)bt.add(a[l].z,-a[l].cnt);
}


int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d %d %d",&b[i].x,&b[i].y,&b[i].z);
    }
    sort(b+1,b+1+n,cmpx);
    int t=0;
    for(int i=1;i<=n;++i){
        t++;
        if(b[i].x!=b[i+1].x||b[i].y!=b[i+1].y||b[i].z!=b[i+1].z){
            a[++top].x=b[i].x;
            a[top].y=b[i].y;
            a[top].z=b[i].z;
            a[top].cnt=t;
            t=0;
        }
    }
	//三维偏序
    cdq(1,top);
    for(int i=1;i<=top;++i)
        ans[a[i].ans+a[i].cnt-1]+=a[i].cnt;
    for(int i=0;i<n;++i)printf("%d\n",ans[i]);
}
