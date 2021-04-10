/*
整体二分与cdq分治有异曲同工之妙。

二分答案的过程中即为对所有询问分治的过程。
*/

//https://www.luogu.com.cn/problem/P1527
/*
N=500,给出N*N的矩阵
1e5次询问，每次询问一个矩形子区间中第k大的值。
*/
#include<bits/stdc++.h>
using namespace std;
const int M=60010;
const int N=510;
int n,m,top,ans[M];

//要分治 先排序
struct nd{
    int x,y,v;
    bool operator < (const nd& t)const {
        return v<t.v;
    }
}a[N*N];
struct Q{
    int x1,y1,x2,y2,k;
    int id;
};
typedef vector<Q> VQ;

struct BIT2{
    int a[N][N];
    int lowbit(int x){return x&-x;}
    void add(int x,int y,int v){
        for(;x<=n;x+=lowbit(x))
            for(int t=y;t<=n;t+=lowbit(t))
                a[x][t]+=v;
    }
    int sum(int x,int y){
        int res=0;
        for(;x>0;x-=lowbit(x))
            for(int t=y;t>0;t-=lowbit(t))
                res+=a[x][t];
        return res;
    }
    int get(int x1,int y1,int x2,int y2){
        int res=sum(x1-1,y1-1)+sum(x2,y2);
        res-=sum(x2,y1-1)+sum(x1-1,y2);
        return res;
    }
}tr;

void cdq(int l,int r,VQ vc){
    if(vc.size()==0)return;
    if(l==r){
        for(auto v:vc)ans[v.id]=a[l].v;
        return;
    }
    int m=l+r>>1;
    for(int i=l;i<=m;++i)
        tr.add(a[i].x,a[i].y,1);
    VQ v1,v2;
    for(auto v:vc){
        int t=tr.get(v.x1,v.y1,v.x2,v.y2);
        if(v.k<=t)v1.push_back(v);
        else{
            v.k-=t;
            v2.push_back(v);
        }
    }
    for(int i=l;i<=m;++i)
        tr.add(a[i].x,a[i].y,-1);

    cdq(l,m,v1);
    cdq(m+1,r,v2);
}

int main(){
    scanf("%d %d" ,&n,&m);
    for(int i=1;i<=n;++i){
        for(int j=1;j<=n;++j){
            a[++top].x=i;
            a[top].y=j;
            scanf("%d",&a[top].v);
        }
    }
    sort(a+1,a+1+top);
    VQ vc;
    Q t;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d %d %d",&t.x1,&t.y1,&t.x2,&t.y2,&t.k);
        t.id=i;
        vc.push_back(t);
    }
    //cdq离线处理
    cdq(1,top,vc);
    for(int i=1;i<=m;++i){
        printf("%d\n",ans[i]);
    }
}
