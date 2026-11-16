// Independent Magma audit of maps.txt, including Hamilton--Loeffler tests.
// Relations: Hamilton--Loeffler, LMS J. Comput. Math. 17 (2014), 206--208.
// The mixed-modulus inverse of 5 is taken modulo the 2-primary factor.
SetSeed(1); SetColumns(0);
function IntegersCSV(st)
    return [StringToInteger(x):x in Split(st,",")];
end function;
function CongruenceTest(S,R,N)
    L:=S^-1*R; B:=S^-1*L^-1*S;
    e:=2^Valuation(N,2); m:=N div e;
    if e eq 1 then
        return (B^2*L^(-InverseMod(2,N)))^3 eq Id(Parent(S));
    elif m eq 1 then
        W:=L^20*B^InverseMod(5,N)*L^-4*B^-1;
        V:=L*B^-1*L;
        return V^-1*W*V eq W^-1 and W^-1*B*W eq B^25
            and (W*B^5*L*B^-1*L)^3 eq V^2;
    end if;
    c:=CRT([0,1],[e,m]); d:=CRT([1,0],[e,m]);
    a:=L^c; b:=B^c; l:=L^d; r:=B^d;
    w:=l^20*r^InverseMod(5,e)*l^-4*r^-1;
    v:=a*b^-1*a; u:=l*r^-1*l;
    return a*r eq r*a and v^4 eq Id(Parent(S)) and v^2 eq (b^-1*a)^3
        and v^2 eq (b^2*a^(-InverseMod(2,m)))^3
        and u^-1*w*u eq w^-1 and w^-1*r*w eq r^25
        and u^2 eq (w*r^5*l*r^-1*l)^3;
end function;

// Validate the test on every previously constructed congruence action first.
for line in Split(Read("congruence_actions.txt"),"\n") do
    if #line eq 0 or line[1] eq "#" then continue; end if;
    f:=Split(line,"|"); rr:=IntegersCSV(f[2]); ss:=IntegersCSV(f[3]); n:=#rr;
    S:=Sym(n)![i+1:i in ss]; R:=Sym(n)![i+1:i in rr];
    assert CongruenceTest(S,R,2*Order(S^-1*R));
end for;
print "PASS: congruence criterion on all 64 original representatives";
total:=0; cong:=0; actual:=0; congactual:=0; projectivecong:=0; overcong:=0;
gcounts:=[0,0]; gc:=[0,0]; counts:=AssociativeArray();
PrintFile("verified-lifts.txt","# map_id|bits|congruence|SL2_index|genus|regular_widths|irregular_widths":Overwrite:=true);
PrintFile("verified-projective.txt","# map_id|congruence|generalized_level":Overwrite:=true);
for line in Split(Read("maps.txt"),"\n") do
    if #line eq 0 or line[1] eq "#" then continue; end if;
    f:=Split(line,"|"); id,mu,e3,g,c,nreg,aut,nraw,nclasses:=Explode([StringToInteger(f[i]):i in [1..9]]);
    rr:=[i+1:i in IntegersCSV(f[10])]; ss:=[i+1:i in IntegersCSV(f[11])]; bitslist:=IntegersCSV(f[12]);
    S:=Sym(mu)!ss; R:=Sym(mu)!rr; P:=sub<Sym(mu)|S,R>;
    assert IsTransitive(P) and S^2 eq Id(P) and R^3 eq Id(P);
    assert forall{i:i in [1..mu]|i^S ne i};
    assert #[i:i in [1..mu]|i^R eq i] eq e3;
    cusps:=Orbits(sub<P|S*R>); widths:=[#o:o in cusps];
    assert #cusps eq c and 1+mu/12-e3/3-c/2 eq g;
    assert nreg eq (mu+2*e3-12)/3;
    A:=Centralizer(Sym(mu),P); assert #A eq aut;
    V:=Orbits(sub<P|R>); Sort(~V,func<a,b|Min(a)-Min(b)>);
    vertex:=[0:i in [1..mu]]; tree:={Integers()|}; starts:=[Min(v):v in V];
    for j in [1..#V] do for i in V[j] do vertex[i]:=j; end for; end for;
    for j in [2..#V] do d:=starts[j]; assert vertex[ss[d]] lt j; tree join:={d,ss[d]}; end for;
    free:=[i:i in [1..mu]|i lt ss[i] and not i in tree];
    assert #free eq 2*g+c-1;
    raw:=[]; voltages:=AssociativeArray();
    for bits in [0..2^#free-1] do
        voltage:=[Integers()|i gt ss[i] select 1 else 0:i in [1..mu]];
        for j in [1..#free] do
            if (bits div 2^(j-1)) mod 2 eq 1 then
                i:=free[j]; voltage[i]:=1-voltage[i]; voltage[ss[i]]:=1-voltage[ss[i]];
            end if;
        end for;
        regular:=[o:o in cusps|(&+[1+voltage[i]:i in o]) mod 2 eq 0];
        if #regular eq nreg then Append(~raw,bits);voltages[bits]:=voltage;end if;
    end for;
    assert #raw eq nraw and nraw eq 2^(2*g)*Binomial(c,nreg);
    // Burnside's lemma gives a second computation of the lift-orbit count.
    fixed:=0; orbitmins:={Integers()|};
    for bits in raw do
        voltage:=voltages[bits]; images:={Integers()|};
        for a in A do
            gauge:=[0:j in [1..#V]];
            for j in [2..#V] do d:=starts[j];gauge[j]:=(gauge[vertex[ss[d]]]+voltage[d^a]+1) mod 2;end for;
            newbits:=&+[Integers()|((voltage[i^a]+gauge[vertex[i]]+gauge[vertex[ss[i]]]) mod 2)*2^(j-1):j in [1..#free]|true where i is free[j]];
            assert newbits in raw;
            Include(~images,newbits);if newbits eq bits then fixed+:=1;end if;
        end for;
        Include(~orbitmins,Min(images));
    end for;
    assert fixed div aut eq nclasses and fixed mod aut eq 0;
    assert Set(bitslist) eq orbitmins;
    N:=2*LCM(widths); pc:=CongruenceTest(S,R,N);
    PrintFile("verified-projective.txt",Sprintf("%o|%o|%o",id,pc,N div 2));
    if pc then projectivecong+:=1;overcong+:=(mu div aut)*nraw;end if;
    actual+:=(mu div aut)*nraw;
    for bits in bitslist do
        voltage:=voltages[bits];
        sl:=Sym(2*mu)![2*(ss[(j div 2)+1]-1)+((j mod 2+voltage[(j div 2)+1]) mod 2)+1:j in [0..2*mu-1]];
        rl:=Sym(2*mu)![2*(rr[(j div 2)+1]-1)+(j mod 2)+1:j in [0..2*mu-1]];
        z:=sl^2;assert Order(sl) eq 4 and rl^3 eq Id(Parent(sl)) and z*rl eq rl*z;
        yes:=CongruenceTest(sl,rl,N);assert not yes or pc;
        reg:=[#o:o in cusps|(&+[1+voltage[i]:i in o]) mod 2 eq 0];
        irr:=[#o:o in cusps|(&+[1+voltage[i]:i in o]) mod 2 eq 1];
        assert mu/6+e3/3-#reg/2 eq 2;
        PrintFile("verified-lifts.txt",Sprintf("%o|%o|%o|%o|%o|%o|%o",id,bits,yes,2*mu,g,Sort(reg),Sort(irr)));
        total+:=1;gcounts[g+1]+:=1;
        key:=<2*mu,g>;if not IsDefined(counts,key) then counts[key]:=[0,0];end if;
        vals:=counts[key];vals[1]+:=1;
        if yes then
            cong+:=1;gc[g+1]+:=1;vals[2]+:=1;
            congactual+:=(2*mu) div #Centralizer(Sym(2*mu),sub<Sym(2*mu)|sl,rl>);
        end if;
        counts[key]:=vals;
    end for;
    if id mod 1000 eq 0 then print id,total,cong;end if;
end for;
print "TOTAL CLASSES",total,"GENUS",gcounts,"ACTUAL",actual;
print "CONGRUENCE CLASSES",cong,"GENUS",gc,"ACTUAL",congactual;
print "PROJECTIVE CONGRUENCE",projectivecong,"ACTUAL LIFTS ABOVE THEM",overcong;
print "NONCONGRUENCE CLASSES",total-cong,"ACTUAL",actual-congactual;
assert total eq 47244 and actual eq 1466010;
assert cong eq 64 and congactual eq 558 and gc eq [32,32];
assert projectivecong eq 44 and overcong eq 1745;
PrintFile("verified-counts.txt","# SL2_index|genus|all_classes|congruence_classes|noncongruence_classes":Overwrite:=true);
for key in Sort(Setseq(Keys(counts))) do
    vals:=counts[key];PrintFile("verified-counts.txt",Sprintf("%o|%o|%o|%o|%o",key[1],key[2],vals[1],vals[2],vals[1]-vals[2]));
end for;
print "PASS: complete census, independent congruence test, Burnside audit";
quit;
