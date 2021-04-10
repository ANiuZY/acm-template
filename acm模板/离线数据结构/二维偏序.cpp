
/*
二维偏序问题：离线问题树状数组，在线问题主席树。

二分套二维偏序问题为在线二维偏序，用主席树二维偏序（或整体二分）可解。
*/

//归并排序求逆序对
#include<bits/stdc++.h>
using namespace std;
const int M=500010;
int a[M],b[M],n;
typedef long long ll;
ll ans;

void Bmerge(int l,int m,int r){
    int x=l,y=m+1,p=l;
    while(p<=r){
        if(y>r||a[x]<=a[y]&&x<=m)b[p++]=a[x++];
        else if(x>m||a[x]>a[y])b[p++]=a[y++],
        ans+=m+1-x;//<-
    }
    while(l<=r)a[l]=b[l],l++;
}

void Bsort(int l,int r){
    if(l==r)return;
    int m=l+r>>1;
    Bsort(l,m);
    Bsort(m+1,r);
    Bmerge(l,m,r);
}

int main(){
    scanf("%d",&n);
    for(int i=1;i<=n;++i)scanf("%d",a+i);
    Bsort(1,n);
    printf("%lld",ans);
}
