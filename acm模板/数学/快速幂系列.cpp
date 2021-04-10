#include<bits/stdc++.h>
using namespace std;
const int M=1000010;
typedef long long ll;
int mod,x0,x1,a,b;
char s[M];

//二进制快速幂
template<class T> T Bpow(T a,ll b){
    T res=1;
    while(b){
        if(b&1)res=(res*a)%mod;
        b>>=1;
        a=(a*a)%mod;
    }
    return res;
}

///十进制快速幂
template<class T> T Spow(T a,char s[]){
    T res=1;
    int p=strlen(s)-1;
    while(p>=0){
        ll k=s[p]-'0';
        res=(res*Bpow(a,k))%mod;
        --p;
        a=Bpow(a,10);
    }
    return res;
}

///矩阵
struct matrix{
    int size;
    ll data[8][8];
    void set(int k){size=k;}
    matrix(){
        size=2;
        memset(data,0,sizeof data);
    }
    matrix(int t){
        size=2;
        memset(data,0,sizeof data);
        for(int i=0;i<size;++i)data[i][i]=t;
    }
    matrix operator * (const matrix& t){
        matrix tmp;
        for(int i=0;i<size;++i)
            for(int j=0;j<size;++j)
                for(int k=0;k<size;++k)
                    tmp.data[i][j]+=data[i][k]*t.data[k][j];
        return tmp;
    }
    matrix operator % (ll mod){
        for(int i=0;i<size;++i)for(int j=0;j<size;++j)data[i][j]%=mod;
        return *this;
    }
    void print(){
        for(int i=0;i<size;++i){
        for(int j=0;j<size;++j)
            cout<<data[j][i]<<' ';
            cout<<endl;
        }
    }
};

/// 牛客多校2019-5-B:
int main(){
    scanf("%d%d%d%d",&x0,&x1,&a,&b);
    scanf("%s%d",s,&mod);
    matrix m;
    m.data[0][0]=a;
    m.data[1][0]=b;
    m.data[0][1]=1;
    matrix ans=Spow(m,s);
    //ans.print();
    printf("%d\n",(ans.data[0][1]*x1%mod+ans.data[1][1]*x0%mod)%mod);
}
