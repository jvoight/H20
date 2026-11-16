SetSeed(1); SetColumns(0);
function MatPerm(N,a)
    return Sym(N^2)![1+(((i mod N)*a[1]+(i div N)*a[3]) mod N)
        +N*(((i mod N)*a[2]+(i div N)*a[4]) mod N):i in [0..N^2-1]];
end function;
function Model(N)
    return sub<Sym(N^2)|MatPerm(N,[0,1,-1,0]),MatPerm(N,[1,1,0,1])>;
end function;
function CSV(p,n)
    st:="";
    for i in [1..n] do st cat:= (i eq 1 select "" else ",") cat IntegerToString(i^p-1); end for;
    return st;
end function;
ChangeDirectory("sources");
load "csg1-lev48.dat";
ChangeDirectory("..");
PrintFile("cp_projective.txt","# CP|Nprime|mu|R|S":Overwrite:=true);
for rec in L do
    m:=rec`index; e3:=#[x:x in rec`c3|x eq 1];
    if m gt 36 or 1 in rec`c2 or m+8*e3 gt 36-12*rec`genus or m+2*e3 lt 12 then continue; end if;
    N:=rec`level; G:=Model(N); d:=MatPerm(N,[-1,0,0,1]);
    J:=sub<G|[MatPerm(N,[Integers()!x:x in a]):a in rec`matgens] cat [G.1^2]>;
    for refl in [0..rec`GLZConj-1] do
        a,P:=CosetAction(G,J^(d^refl));
        s:=a(G.1); r:=a(G.1*G.2);
        assert Degree(P) eq m and Order(s) eq 2 and Order(r) eq 3;
        PrintFile("cp_projective.txt",Sprintf("%o:%o|%o|%o|%o|%o",rec`name,refl,N,m,CSV(r,m),CSV(s,m)));
    end for;
end for;
load "../congruence-enumeration/representatives.m";
PrintFile("congruence_actions.txt","# old_id|R|S":Overwrite:=true);
for id in [1..#CongruencePG2] do
    rec:=CongruencePG2[id]; N:=rec[3]; n:=2*rec[5]; G:=Model(N);
    H:=sub<G|[MatPerm(N,a):a in rec[11]]>;
    a,P:=CosetAction(G,H); s:=a(G.1); r:=a(G.1*G.2);
    assert Degree(P) eq n and Order(s) eq 4 and Order(r) eq 3;
    PrintFile("congruence_actions.txt",Sprintf("%o|%o|%o",id,CSV(r,n),CSV(s,n)));
end for;
print "PASS: exported projective congruence classes and 64 signed actions";
quit;
