// Run after enumerate_lifts.m and fetch_lmfdb.py.
SetSeed(1);
SetColumns(0);
load "lift_representatives.m";
load "lmfdb_candidates.m";
function MatPerm(N,a)
    return Sym(N^2)![1+(((i mod N)*a[1]+(i div N)*a[3]) mod N)
        +N*(((i mod N)*a[2]+(i div N)*a[4]) mod N):i in [0..N^2-1]];
end function;
function Entries(N,p)
    u:=2^p-1; v:=(N+1)^p-1;
    return [u mod N,u div N,v mod N,v div N];
end function;
function Models(N)
    S:=sub<Sym(N^2)|MatPerm(N,[0,1,-1,0]),MatPerm(N,[1,1,0,1])>;
    E:=sub<Sym(N^2)|S,[MatPerm(N,[u,0,0,1]):u in [1..N-1]|GCD(u,N) eq 1]>;
    return S,E;
end function;
models:=AssociativeArray();
dbgroups:=AssociativeArray();
for r in LMFDB do
    N:=r[2];
    if not IsDefined(models,N) then S,E:=Models(N); models[N]:=<S,E>; end if;
    S:=models[N][1]; E:=models[N][2];
    G:=sub<E|[MatPerm(N,a):a in r[6]]>;
    assert Index(E,G) eq r[3];
    assert #G div #(G meet S) eq EulerPhi(N);
    assert not MatPerm(N,[-1,0,0,-1]) in G;
    dbgroups[r[1]]:=G;
end for;
PrintFile("matches.txt","# id|modulus|Hlevel|lift_number|Glevel|fine_label|coarse_label|conjugator":Overwrite:=true);
for lift in RationalLifts do
    id,N,hlevel,k,gens:=Explode(lift);
    S:=models[N][1]; E:=models[N][2];
    G:=sub<E|[MatPerm(N,a):a in gens]>;
    level:=N;
    if N ne hlevel then
        S0:=models[hlevel][1]; E0:=models[hlevel][2];
        G0:=sub<E0|[MatPerm(hlevel,a):a in gens]>;
        if #G eq #G0*(#E div #E0) then
            G:=G0; E:=E0; S:=S0; level:=hlevel;
        end if;
    end if;
    idx:=Index(E,G);
    matches:=[]; conj:=[];
    for r in LMFDB do
        if r[2] ne level or r[3] ne idx then continue; end if;
        yes,x:=IsConjugate(E,G,dbgroups[r[1]]);
        if yes then
            assert G^x eq dbgroups[r[1]];
            Append(~matches,r); conj:=Entries(level,x);
        end if;
    end for;
    assert #matches eq 1;
    r:=matches[1];
    line:=Sprintf("%o|%o|%o|%o|%o|%o|%o|%o",id,N,hlevel,k,level,r[1],r[5],conj);
    PrintFile("matches.txt",line);
    print line;
end for;
print "PASS: every lift matched to exactly one LMFDB fine label by GL2 conjugacy";
quit;
