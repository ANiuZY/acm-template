//https://www.luogu.com.cn/problem/P4396
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
int n,m,a[M],pre[M],lst[M],ans1[M],ans2[M];

struct bit{
    int a[M];
    int lowbit(int x){return x&-x;};
    void add(int x,int val){
        for(;x<M;x+=lowbit(x))
            a[x]+=val;
    }
    int get(int x){
        int res=0;
        for(;x>0;x-=lowbit(x))
            res+=a[x];
        return res;
    }
}tr;

struct nd{
    int x,y,z,f,id;
    bool operator < (const nd& t)const{
        return y<t.y;
    }
}b[M*5];
int cnt;
bool cmp(const nd& t1,const nd& t2){
    if(t1.x==t2.x)return t1<t2;
    return t1.x<t2.x;
}

void cdq(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    cdq(l,m);
    cdq(m+1,r);

    int x=l,y=m+1,t=0;
    while(y<=r){
        if(x<=m&&b[x].y<=b[y].y){
            if(b[x].f==0)tr.add(b[x].z+1,1),t++;
            x++;
        }else{
            if(b[y].f!=0){
                ans1[b[y].id]+=t*b[y].f;
                ans2[b[y].id]+=tr.get(b[y].z+1)*b[y].f;
            }
            y++;
        }
    }
    for(int i=l;i<x;++i)
        if(b[i].f==0)tr.add(b[i].z+1,-1);
    stable_sort(b+l,b+r+1);
    //归并排序中stable_sort比sort快。
}

int main(){
    scanf("%d %d" ,&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        pre[i]=lst[a[i]];
        lst[a[i]]=i;
        b[++cnt]=(nd){i,a[i],pre[i],0,0};
    }
    int l,r,A,B;
    for(int i=1;i<=m;++i){
        scanf("%d %d %d %d",&l,&r,&A,&B);//询问拆解
        b[++cnt]=(nd){r,B,l-1,1,i};
        b[++cnt]=(nd){l-1,A-1,l-1,1,i};
        b[++cnt]=(nd){r,A-1,l-1,-1,i};
        b[++cnt]=(nd){l-1,B,l-1,-1,i};
    }
    stable_sort(b+1,b+cnt+1,cmp);
    cdq(1,cnt);
    for(int i=1;i<=m;++i)printf("%d %d\n",ans1[i],ans2[i]);
}
