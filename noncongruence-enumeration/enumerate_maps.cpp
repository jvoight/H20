// Exhaustive rooted oriented maps with vertices of degree 1 or 3.
// S pairs darts and R rotates vertices.  S has no fixed points.
// New vertices receive canonical breadth-first labels, so every rooted map
// occurs once.  Rerooting removes conjugacy duplicates and computes Aut.
// Signed actions give all lifts to SL2(Z) avoiding -I, modulo Aut.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
constexpr int MAX=36;
int m, want1, want3, s[MAX], r[MAX], vert[MAX], starts[MAX];
uint64_t leaves, rooted[2], maps[2], liftclasses[2], actual[2];
ofstream records("maps.txt"), summary("summary.csv");
int recordid=0;
uint64_t choose(int n,int k){uint64_t z=1;for(int i=1;i<=k;++i)z=z*(n-i+1)/i;return z;}

// Return the canonical labeling based at root.  p maps new labels to old.
void reroot(int root,array<int,MAX>& ns,array<int,MAX>& nr,array<int,MAX>& p){
    int label[MAX];fill(label,label+m,-1);int used=0;
    auto addvertex=[&](int a){int b=a;do{label[b]=used;p[used++]=b;b=r[b];}while(b!=a);};
    addvertex(root);
    for(int i=0;i<used;++i){int a=p[i];if(label[s[a]]<0)addvertex(s[a]);ns[i]=label[s[a]];nr[i]=label[r[a]];}
    assert(used==m);
}
void finish(int nv){
    ++leaves;
    bool seen[MAX]={};vector<vector<int>> cusps;
    for(int i=0;i<m;++i)if(!seen[i]){vector<int> c;int a=i;do{seen[a]=true;c.push_back(a);a=r[s[a]];}while(a!=i);cusps.push_back(c);}
    int c=cusps.size(), numer=2-nv+m/2-c;
    assert(numer>=0 && numer%2==0);int g=numer/2;
    int reg=(m+2*want1-12)/3;
    if(g>1 || reg<0 || reg>c || reg%2)return;
    ++rooted[g];
    // For maps with leaves, a leaf-rooted labeling is lexicographically first.
    if(want1 && r[0]!=0)return;
    vector<array<int,MAX>> autos;
    array<int,MAX> ns,nr,p;
    for(int root=0;root<m;++root){
        if(want1 && r[root]!=root)continue;
        reroot(root,ns,nr,p);
        int cmp=0;
        for(int i=0;i<m && !cmp;++i)cmp=(nr[i]>r[i])-(nr[i]<r[i]);
        for(int i=0;i<m && !cmp;++i)cmp=(ns[i]>s[i])-(ns[i]<s[i]);
        if(cmp<0)return;
        if(cmp==0)autos.push_back(p);
    }
    assert(!autos.empty());int aut=autos.size();assert(m%aut==0);
    ++maps[g];
    // R-voltages are gauged to zero.  Gauge the discovery-tree edges to
    // S-voltage 0 at the earlier dart and 1 at the later dart.
    bool tree[MAX]={};
    for(int v=1;v<nv;++v){int d=starts[v];assert(vert[s[d]]<v);tree[d]=tree[s[d]]=true;}
    vector<int> freeedges;
    for(int i=0;i<m;++i)if(i<s[i]&&!tree[i])freeedges.push_back(i);
    int b=freeedges.size();assert(b==2*g+c-1);
    vector<int> good, representatives;
    for(int bits=0;bits<(1<<b);++bits){
        int voltage[MAX];for(int i=0;i<m;++i)voltage[i]=(i>s[i]);
        for(int j=0;j<b;++j)if((bits>>j)&1){voltage[freeedges[j]]^=1;voltage[s[freeedges[j]]]^=1;}
        int regular=0;
        for(auto& cusp:cusps){int parity=0;for(int i:cusp)parity^=1^voltage[i];regular+=(parity==0);}
        if(regular!=reg)continue;
        good.push_back(bits);int smallest=bits;
        for(auto& a:autos){
            int gauge[MAX]={};
            for(int v=1;v<nv;++v){int d=starts[v];gauge[v]=gauge[vert[s[d]]]^voltage[a[d]]^1;}
            int image=0;
            for(int j=0;j<b;++j){int i=freeedges[j];image|=(voltage[a[i]]^gauge[vert[i]]^gauge[vert[s[i]]])<<j;}
            smallest=min(smallest,image);
        }
        if(bits==smallest)representatives.push_back(bits);
    }
    assert(good.size()==(1u<<(2*g))*choose(c,reg));
    liftclasses[g]+=representatives.size();actual[g]+=(m/aut)*good.size();
    records<<++recordid<<'|'<<m<<'|'<<want1<<'|'<<g<<'|'<<c<<'|'<<reg<<'|'<<aut<<'|'<<good.size()<<'|'<<representatives.size();
    for(auto* ar:{r,s}){records<<'|';for(int i=0;i<m;++i){if(i)records<<',';records<<ar[i];}}
    records<<'|';for(size_t i=0;i<representatives.size();++i){if(i)records<<',';records<<representatives[i];}
    records<<'\n';
}
void generate(int used,int n1,int n3,int nv){
    int x=0;while(x<used&&s[x]>=0)++x;
    if(x==used){if(used==m)finish(nv);return;}
    for(int y=x+1;y<used;++y)if(s[y]<0){s[x]=y;s[y]=x;generate(used,n1,n3,nv);s[x]=s[y]=-1;}
    if(n1<want1){int y=used;s[x]=y;s[y]=x;r[y]=y;vert[y]=nv;starts[nv]=y;generate(used+1,n1+1,n3,nv+1);s[x]=s[y]=-1;}
    if(n3<want3){int y=used;s[x]=y;s[y]=x;s[y+1]=s[y+2]=-1;r[y]=y+1;r[y+1]=y+2;r[y+2]=y;vert[y]=vert[y+1]=vert[y+2]=nv;starts[nv]=y;generate(used+3,n1,n3+1,nv+1);s[x]=s[y]=s[y+1]=s[y+2]=-1;}
}
int main(){
    summary<<"mu,e3,genus,projective_actual,projective_classes,lift_actual,lift_classes\n";
    records<<"# id|mu|e3|genus|cusps|regular|automorphisms|raw_lifts|lift_classes|R|S|representative_bits\n";
    for(m=2;m<=36;m+=2)for(want1=0;want1<=m;++want1){
        if((m-want1)%3)continue;want3=(m-want1)/3;
        int regnum=m+2*want1-12;
        if(regnum<0||regnum%3||(regnum/3)%2||m+8*want1>36)continue;
        leaves=0;fill(rooted,rooted+2,0);fill(maps,maps+2,0);fill(liftclasses,liftclasses+2,0);fill(actual,actual+2,0);
        fill(s,s+MAX,-1);starts[0]=0;
        if(want1){r[0]=0;vert[0]=0;generate(1,1,0,1);}
        if(want3){r[0]=1;r[1]=2;r[2]=0;vert[0]=vert[1]=vert[2]=0;generate(3,0,1,1);}
        for(int g=0;g<2;++g)if(maps[g]){
            summary<<m<<','<<want1<<','<<g<<','<<rooted[g]<<','<<maps[g]<<','<<actual[g]<<','<<liftclasses[g]<<'\n';
            cerr<<"mu="<<m<<" e3="<<want1<<" g="<<g<<" projective="<<maps[g]<<" lift_classes="<<liftclasses[g]<<" actual="<<actual[g]<<"\n";
        }
        summary.flush();records.flush();cerr<<"all rooted maps for signature: "<<leaves<<"\n";
    }
    cerr<<"PASS: exhaustive enumeration and lift-count identities\n";
}
