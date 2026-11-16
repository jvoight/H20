// Run here with MAGMA_ARCH=arm64 magma -b enumerate_lifts.m
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
function Models(N)
    S:=sub<Sym(N^2)|MatPerm(N,[0,1,-1,0]),MatPerm(N,[1,1,0,1])>;
    E:=sub<Sym(N^2)|S,[MatPerm(N,[u,0,0,1]):u in [1..N-1]|GCD(u,N) eq 1]>;
    assert Index(E,S) eq EulerPhi(N);
    return S,E;
end function;
PrintFile("higher-lifts.txt","# id|CP|N|exact_level|mode|normalizer_quotient_order|det_normalizer_order|lift_count":Overwrite:=true);
PrintFile("higher_lift_representatives.m","// <id, modulus, exact level of H, lift number, generators>\nRationalLifts := [* *];":Overwrite:=true);
for id in [2,6,7] do
    r:=CongruencePG2[id];
    for N in [2*r[4],4*r[4]] do
        S,E:=Models(N);
        S0,E0:=Models(r[4]);
        red:=hom<S->S0|S0.1,S0.2>;
        H0:=sub<S0|[MatPerm(r[4],a):a in r[11]]>;
        H:=H0@@red;
        assert Index(S,H) eq 2*r[5];
        assert not MatPerm(N,[-1,0,0,-1]) in H;
        A:=Normalizer(E,H);
        Q,q:=quo<A|H>;
        K:=q(A meet S);
        d:=#Q div #K;
        good:=[];
        if d eq EulerPhi(N) then
            for c in Subgroups(Q:OrderEqual:=EulerPhi(N)) do
                C:=c`subgroup;
                if #(C meet K) ne 1 then continue; end if;
                G:=C@@q;
                assert G meet S eq H and #G eq #H*EulerPhi(N);
                assert Index(E,G) eq 2*r[5];
                assert forall{P:P in good|not IsConjugate(E,G,P)};
                Append(~good,G);
                gens:=[Entries(N,G.i):i in [1..Ngens(G)]];
                PrintFile("higher_lift_representatives.m",Sprintf("Append(~RationalLifts,<%o,%o,%o,%o,%o>);",id,N,r[4],#good,gens));
            end for;
        end if;
        mode:=N eq r[4] select "exact" else "proof";
        line:=Sprintf("%o|%o|%o|%o|%o|%o|%o|%o",id,r[1],N,r[4],mode,#Q,d,#good);
        print line;
        PrintFile("higher-lifts.txt",line);
    end for;
end for;
print "PASS: exceptional classes checked at doubled and quadrupled levels";
quit;
