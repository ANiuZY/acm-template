#include<bits/stdc++.h>
#define M 200010
#define mod 998244353
typedef long long ll;
using namespace std;

//o(n)预处理，o(1)查询
ll f[M],inv[M];
ll qpow(ll a,ll b){
    ll res=1;
    while(b){
        if(b&1)res=res*a%mod;
        a=a*a%mod;b>>=1;
    }
    return res;
}
ll C(int n,int m){
    if(n<0||m<0||m>n)return 0;
    if(m==0||m==n)return 1;
    return f[n]*inv[n-m]%mod*inv[m]%mod;
}
void CInit(){
    f[1]=1;
    for(int i=2;i<M;i++)f[i]=(f[i-1]*i)%mod;
    inv[M-1]=qpow(f[M-1],mod-2);
    for(int i=M-2;i>=1;i--)inv[i]=(inv[i+1]*(i+1))%mod;
}

//----------------

//递推预处理 o(n^2)预处理，o(1)查询
int C[M][M];
void init(){
    C[0][0]=1;
    for(int i=1;i<M;++i){
        C[i][0]=1;
        for(int j=1;j<=i;++j)
            C[i][j]=(C[i-1][j]+C[i-1][j-1])%mod;
    }
}
