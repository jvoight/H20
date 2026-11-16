// Run from this directory: magma -b enumerate_congruence.m
SetSeed(1);
SetColumns(0);
ChangeDirectory("sources");
load "csg1-lev48.dat";
ChangeDirectory("..");
CP:=[r:r in L|r`index le 36 and not 1 in r`c2];

// Faithful row-vector action on (Z/mZ)^2. The identity vector is retained.
function MatPerm(m, a)
    return Sym(m^2)![1+(((i mod m)*a[1]+(i div m)*a[3]) mod m)
        +m*(((i mod m)*a[2]+(i div m)*a[4]) mod m): i in [0..m^2-1]];
end function;
function Model(m)
    s:=MatPerm(m,[0,1,-1,0]); t:=MatPerm(m,[1,1,0,1]);
    return sub<Sym(m^2)|s,t>;
end function;
function Entries(m,p)
    u:=(2^p)-1; v:=((m+1)^p)-1;
    return [u mod m,u div m,v mod m,v div m];
end function;
function Word(G,w)
    return &*([G!1] cat [G.Abs(i)^Sign(i):i in w]);
end function;

PrintFile("results.txt","# CP reflection level mu genus e3 regular_widths irregular_widths orbit_size generators_mod_2N":Overwrite:=true);
PrintFile("audit.txt","# CP tested_lifts qualifying_lifts SL_classes":Overwrite:=true);
PrintFile("representatives.m","// <CP label, projective level, modulus, level, mu, genus, e3, regular widths, irregular widths, conjugacy orbit size, matrix generators>\nCongruencePG2 := [* *];":Overwrite:=true);
count:=0; literal:=0; expected_literal:=0; glcount:=0;
for r in CP do
    e3:=#[x:x in r`c3|x eq 1]; mu:=r`index;
    cr:=(mu+2*e3-12)/3;
    if cr lt 0 or cr gt #r`cusps or Denominator(cr) ne 1 then continue; end if;
    N:=r`level; M:=2*N;
    G:=Model(M); B:=Model(N);
    red:=hom<G->B|B.1,B.2>;
    J:=sub<B|[Word(B,w):w in r`generators] cat [B.1^2]>;
    assert Index(B,J) eq mu;
    // Verify the alternate matrix encoding in the published data.
    assert J eq sub<B|[MatPerm(N,[Integers()!x:x in a]):a in r`matgens] cat [B.1^2]>;
    H:=J@@red;
    assert Index(G,H) eq mu;
    z:=G.1^2; t:=G.2; s:=G.1;
    assert z in H and Order(z) eq 2;
    act,Q:=CosetAction(G,H);
    assert #[i:i in [1..mu]|i^act(s) eq i] eq 0;
    assert #[i:i in [1..mu]|i^act(s*t) eq i] eq e3;
    cyc:=Orbits(sub<Q|act(t)>);
    assert Sort([#o:o in cyc]) eq Sort(r`cusps);
    // A transversal records positive unipotent cusp generators.
    trans:=Transversal(G,H);
    par:=[]; widths:=[];
    // Use double cosets via the action, with transversal representatives
    // verified directly by membership (independent of action numbering).
    todo:={x:x in trans};
    while #todo gt 0 do
        x:=Rep(todo); w:=1;
        while not x*t^w*x^-1 in H do w+:=1; end while;
        Append(~par,x*t^w*x^-1); Append(~widths,w);
        todo:={y:y in todo|not exists{k:k in [0..w-1]|y*(x*t^k)^-1 in H}};
    end while;
    assert Sort(widths) eq Sort(r`cusps);
    assert LCM(widths) eq N;
    assert 1+mu/12-e3/3-#widths/2 eq r`genus;
    assert Index(G,Normalizer(G,H)) eq r`length;
    A,a:=AbelianQuotient(H); C:=AbelianGroup([2]);
    good:=[]; tested:=0;
    for bits in CartesianPower([0,1],Ngens(A)) do
        if exists{i:i in [1..Ngens(A)]|bits[i] eq 1 and IsOdd(Order(A.i))} then continue; end if;
        chi:=hom<A->C|[bits[i]*C.1:i in [1..Ngens(A)]]>;
        if chi(a(z)) eq C!0 then continue; end if;
        tested+:=1;
        reg:=[widths[i]:i in [1..#par]|chi(a(par[i])) eq C!0];
        irr:=[widths[i]:i in [1..#par]|chi(a(par[i])) ne C!0];
        if #reg ne cr then continue; end if;
        K:=Kernel(a*chi);
        assert not z in K and Index(G,K) eq 2*mu;
        Append(~good,<K,Sort(reg),Sort(irr)>);
    end for;
    expected_literal+:=r`GLZConj*r`length*#good;
    reps:=[];
    for v in good do
        if not exists{u:u in reps|IsConjugate(G,u[1],v[1])} then Append(~reps,v); end if;
    end for;
    PrintFile("audit.txt",Sprintf("%o %o %o %o",r`name,tested,#good,#reps));
    printf "%o: %o congruence lifts; %o with pg=2; %o SL classes\n",r`name,tested,#good,#reps;
    d:=MatPerm(M,[-1,0,0,1]);
    E:=sub<Sym(M^2)|G,d>;
    greps:=[];
    for v in reps do
        if not exists{u:u in greps|IsConjugate(E,u[1],v[1])} then Append(~greps,v); end if;
    end for;
    glcount+:=#greps;
    for refl in [0..r`GLZConj-1] do
      for v in reps do
        K:=v[1]^(d^refl);
        level:=Kernel(red) subset K select N else M;
        orb:=Index(G,Normalizer(G,K));
        // Independently count signed cusps on the double coset cover G/K.
        ak,Qk:=CosetAction(G,K);
        cycles:=Orbits(sub<Qk|ak(t)>);
        reg2:=[]; irr2:=[];
        for o in cycles do
            if (Rep(o)^ak(z)) in o then
                assert IsEven(#o); Append(~irr2,#o div 2);
            else
                Append(~reg2,#o);
            end if;
        end for;
        assert Sort(reg2) eq Sort(v[2] cat v[2]);
        assert Sort(irr2) eq v[3];
        assert mu/6+e3/3-#reg2/4 eq 2;
        assert mu+6*#irr2+8*e3 eq 12*(3-r`genus);
        // Check minimal congruence level, not merely a defining modulus.
        for ell in Divisors(M) do
            if ell eq 1 then continue; end if;
            D:=Model(ell); rr:=hom<G->D|D.1,D.2>;
            assert (Kernel(rr) subset K) eq (ell mod level eq 0);
        end for;
        gg:=[Entries(M,k):k in OrderedGenerators(K)];
        assert K eq sub<G|[MatPerm(M,x):x in gg]>;
        PrintFile("representatives.m",Sprintf("Append(~CongruencePG2, <\"%o\",%o,%o,%o,%o,%o,%o,%o,%o,%o,%o>);",r`name,N,M,level,mu,r`genus,e3,v[2],v[3],orb,gg));
        count+:=1; literal+:=orb;
        PrintFile("results.txt",Sprintf("%o|%o|%o|%o|%o|%o|%o|%o|%o|%o",r`name,refl,level,mu,r`genus,e3,v[2],v[3],orb,[Entries(M,k):k in OrderedGenerators(K)]));
      end for;
    end for;
end for;
printf "TOTAL SL conjugacy classes: %o; actual subgroups: %o\n",count,literal;
assert literal eq expected_literal;
printf "GL(2,Z) conjugacy classes: %o\n",glcount;
assert count eq 64 and literal eq 558;
print "PASS: exhaustive congruence lift enumeration and independent cusp/level checks";
quit;
