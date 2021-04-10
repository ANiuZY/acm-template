//求最大异或值
#include<bits/stdc++.h>
using namespace std;
const int M=1000010;
const int N=25;
int n,m,a[M],sum[M];
int rt[M],tr[M<<4][2],num[M<<4],ed[M<<4],cnt;

int clone(int o){
    ++cnt;
    tr[cnt][0]=tr[o][0];
    tr[cnt][1]=tr[o][1];
    num[cnt]=num[o];
    return cnt;
}
void add(int& o,int x){//添加一个数x
    o=clone(o);
    int now=o;
    num[now]++;
    for(int i=N-1;i>=0;--i){
        int t=(x>>i)&1;
        tr[now][t]=clone(tr[now][t]);
        now=tr[now][t];
        num[now]++;
    }
    ed[now]=x;
}
int getMax(int lo,int ro,int ex){//lo-ro之间异或ex后的最大的数
    for(int i=N-1;i>=0;--i){
        int t=(ex>>i)&1;
        int x=num[tr[ro][t^1]]-num[tr[lo][t^1]];
        if(x>0){
            ro=tr[ro][t^1];
            lo=tr[lo][t^1];
        }else{
            ro=tr[ro][t^0];
            lo=tr[lo][t^0];
        }
    }
    return ed[ro];
}
ll getK(int o,ll ex,int k){//字典树中异或ex后k大数
    for(int i=N-1;i>=0;--i){
        int t=(ex>>i)&1;
        int x=num[tr[o][1^t]];
        if(x>=k)o=tr[o][1^t];
        else o=tr[o][0^t],k-=x;
    }
    return ed[o]^ex;
}

int main(){
    scanf("%d %d",&n,&m);
    for(int i=1;i<=n;++i){
        scanf("%d",a+i);
        sum[i]=sum[i-1]^a[i];
        rt[i]=rt[i-1];
        add(rt[i],sum[i]);
    }
    char s[5];
    int x,l,r;
    for(int i=1;i<=m;++i){
        scanf("%s",s);
        if(s[0]=='A'){
            scanf("%d",&x);
            ++n;
            sum[n]=sum[n-1]^x;
            rt[n]=rt[n-1];
            add(rt[n],sum[n]);
        }else{
            scanf("%d %d %d",&l,&r,&x);
            printf("%d\n",getMax(rt[max(0,l-2)],rt[r-1],x^sum[n])^x^sum[n]);
        }
    }
}
