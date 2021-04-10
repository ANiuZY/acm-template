//不考虑交叉点;李超线段树;
#include<bits/stdc++.h>
using namespace std;
const int mod=39989; //x轴范围
const int mod1=1e9;  //y轴范围
int n,ans,cnt;
double h;

struct nd{
    int id;
    double y1,y2;
    nd(){}
    nd(double x1,double x2,int x3){y1=x1;y2=x2;id=x3;}
    double get(int x,int l,int r){
        if(l==r)return max(y1,y2);
        return y1+(y2-y1)*(x-l)/(r-l);
    }
}a[mod<<2];

#define ls o<<1
#define rs o<<1|1
nd k;
int L,R,p;
void add(int o,int l,int r){ // 添加一条[L,R]内的线段k
    if(l>R||r<L)return;
    if(l>=L&&r<=R){
        double y1=k.get(l,L,R),y2=k.get(r,L,R);
        if(y1<a[o].y1&&y2<a[o].y2)return;
        if(y1>=a[o].y1&&y2>=a[o].y2){
            a[o]=nd(y1,y2,k.id);
            return;
        }
    }
    int m=l+r>>1;
    add(ls,l,m);
    add(rs,m+1,r);
}
void get(int o,int l,int r){ // x=p处的最大值
    double t=a[o].get(p,l,r);
    if(t>h)ans=a[o].id,h=t;
    if(l==r)return;
    int m=l+r>>1;
    if(p<=m)get(ls,l,m);
    else get(rs,m+1,r);
}

int main(){
    scanf("%d",&n);
    int op,y1,y2;
    while(n--){
        scanf("%d",&op);
        if(op){
            scanf("%d %d %d %d",&L,&y1,&R,&y2);//强制在线
            k.id=++cnt;
            L=(L+ans-1)%mod+1;
            R=(R+ans-1)%mod+1;
            k.y1=(y1+ans-1)%mod1+1;
            k.y2=(y2+ans-1)%mod1+1;
            if(L>R)swap(L,R),swap(k.y1,k.y2);
            add(1,1,mod);
        }else{
            scanf("%d",&p);
            p=(p+ans-1)%mod+1;
            ans=0;h=0;
            get(1,1,mod);
            printf("%d\n",ans);
        }
    }
}
