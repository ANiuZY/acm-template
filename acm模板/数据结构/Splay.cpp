//单点插入、删除 前驱 后继 第k大 求排名
#include<bits/stdc++.h>
using namespace std;
const int M=100010;
const int inf=1<<30;

struct Splay{
    int x[M]={-inf},fa[M],ls[M],rs[M],sum[M],num[M];
    int& root=rs[0];
    int cnt;
    void pu(int o){
        sum[o]=sum[ls[o]]+sum[rs[o]]+num[o];
    }
    int getdir(int o){
        return ls[fa[o]]==o?0:1;
    }
    void connect(int o,int f,int dir){
        fa[o]=f;
        if(dir)rs[f]=o;
        else ls[f]=o;
    }
    void rotate(int o){
        int f=fa[o];
        int ff=fa[f];
        int dirf=getdir(f);
        int dirx=getdir(o);
        int t=(dirx^1)?rs[o]:ls[o];
        connect(t,f,dirx);connect(f,o,dirx^1);connect(o,ff,dirf);
        pu(f);pu(o);
    }
    void splay(int o,int to){
        to=fa[to];
        while(fa[o]!=to){
            int f=fa[o];
            if(fa[f]==to)rotate(o);
            else if(getdir(o)==getdir(f)){
                rotate(f);
                rotate(o);
            }else{
                rotate(o);
                rotate(o);
            }
        }
    }
    int newpoint(int v,int f){
        x[++cnt]=v;
        fa[cnt]=f;
        sum[cnt]=num[cnt]=1;
        return cnt;
    }
    int find(int v){
        int now=root;
        while(1){
            if(x[now]==v){
                splay(now,root);
                return now;
            }
            now=v<x[now]?ls[now]:rs[now];
            if(now==0)return 0;
        }
    }
    ///---
    void add(int v){
        int now=root;
        while(1){
            if(now)sum[now]++;
            if(v==x[now]){
                num[now]++;
                break;
            }
            int& t=v<x[now]?ls[now]:rs[now];
            if(t==0){
                now=t=newpoint(v,now);
                break;
            }
            now=t;
        }
        splay(now,root);
    }
    void remove(int v){
        int o=find(v);
        if(o==0)return;
        if(num[o]>1){
            num[o]--;
            sum[o]--;
            return;
        }
        if(ls[o]==0){
            root=rs[o];
            fa[root]=0;
        }else{
            int l=ls[o];
            while(rs[l])l=rs[l];
            splay(l,ls[o]);
            connect(rs[o],l,1);
            connect(l,0,1);
            pu(l);
        }
    }
    int getRank(int v){
        int res=0,o=root;
        while(1){
            if(o==0)break;
            if(x[o]==v){
                res+=sum[ls[o]]+1;
                break;
            }
            if(v<x[o])o=ls[o];
            else {
                res+=sum[ls[o]]+num[o];
                o=rs[o];
            }
        }
        if(o)splay(o,root);
        return res;
    }
    int rank(int p){
        if(p<=0||p>sum[root])return -inf;
        int o=root;
        while(1){
            int m=sum[o]-sum[rs[o]];
            if(p>sum[ls[o]] && p<=m)break;
            if(p<m)o=ls[o];
            else {
                p-=m;
                o=rs[o];
            }
        }
        splay(o,root);
        return x[o];
    }
    int upper(int v){
        int o=root,res=inf;
        while(o){
            if(x[o]>v && x[o]<res)res=x[o];
            if(v<x[o])o=ls[o];
            else o=rs[o];
        }
        return res;
    }
    int lower(int v){
        int o=root,res=-inf;
        while(o){
            if(x[o]<v && x[o]>res)res=x[o];
            if(v>x[o])o=rs[o];
            else o=ls[o];
        }
        return res;
    }
};
