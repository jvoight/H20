SetSeed(1);
SetColumns(0);
load "../congruence-enumeration/representatives.m";
function MatPerm(N,a)
    return Sym(N^2)![1+(((i mod N)*a[1]+(i div N)*a[3]) mod N)
        +N*(((i mod N)*a[2]+(i div N)*a[4]) mod N):i in [0..N^2-1]];
end function;
function Entries(N,p)
    u:=2^p-1; v:=(N+1)^p-1;
    return [u mod N,u div N,v mod N,v div N];
end function;
for id in [2,6,7] do
    r:=CongruencePG2[id]; N:=r[4];
    S:=sub<Sym(N^2)|MatPerm(N,[0,1,-1,0]),MatPerm(N,[1,1,0,1])>;
    E:=sub<Sym(N^2)|S,[MatPerm(N,[u,0,0,1]):u in [1..N-1]|GCD(u,N) eq 1]>;
    H:=sub<S|[MatPerm(N,a):a in r[11]]>;
    A:=Normalizer(E,H); Q,q:=quo<A|H>; K:=q(A meet S);
    print "ID",id,"N",N,"H",GroupName(H),"A",GroupName(A),"Q",GroupName(Q),"K",GroupName(K);
    U,u:=UnitGroup(Integers(N));
    dets:=[];
    for x in Q do
        a:=Entries(N,x@@q); d:=(a[1]*a[4]-a[2]*a[3]) mod N;
        Append(~dets,<x,d>);
    end for;
    for d in [d:d in [1..N-1]|GCD(d,N) eq 1] do
        print d,Sort([Order(v[1]):v in dets|v[2] eq d]);
    end for;
    // An involution with determinant -1 is necessary for real descent.
    print "involutions det -1",#[v:v in dets|v[2] eq N-1 and Order(v[1]) eq 2];
    abfull:=[c`subgroup:c in Subgroups(Q)|IsAbelian(c`subgroup) and #(c`subgroup) div #(c`subgroup meet K) eq EulerPhi(N)];
    print "abelian det-surjective subgroups",[GroupName(C):C in abfull];
    realcounts:=[#[v:v in dets|v[1] in C and v[2] eq N-1 and Order(v[1]) eq 2]:C in abfull];
    print "det -1 involutions in these abelian subgroups",realcounts;
    assert forall{n:n in realcounts|n eq 0};
end for;
print "PASS";
quit;
