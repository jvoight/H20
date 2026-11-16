// Run from verification/: magma -b modular_dimensions.m
// Reconstructed dimension checks cited in rewrite-v11.tex.
// Full cuspidal modular symbols (sign 0) have twice the cusp-form dimension.

// Count projective Gamma_1(N) cosets using primitive bottom rows modulo N,
// identifying v with -v. Right translation by T gives the cusp-width orbits.
function Gamma1CuspWidths(N)
    assert N ge 5;
    function CanonicalPair(c,d)
        v := <c mod N,d mod N>;
        w := <-c mod N,-d mod N>;
        return v lt w select v else w;
    end function;
    remaining := {CanonicalPair(c,d) : c,d in [0..N-1] | GCD(GCD(c,d),N) eq 1};
    projective_index := #remaining;
    widths := [];
    while #remaining gt 0 do
        v := Rep(remaining); w := v; width := 0;
        repeat
            assert w in remaining;
            Exclude(~remaining,w);
            w := CanonicalPair(w[1],w[1]+w[2]);
            width +:= 1;
        until w eq v;
        Append(~widths,width);
    end while;
    assert &+widths eq projective_index;
    return Sort(widths),projective_index;
end function;

procedure VerifyModularDimensions(LevelBound)
    assert LevelBound ge 18;
    SetSeed(1);
    SetColumns(0);
    v1,v2,v3 := GetVersion();
    printf "Modular dimensions, Magma %o.%o-%o; levels 1..%o\n",v1,v2,v3,LevelBound;
    small := [0,0,0,0,0,0,1,1,2,2,5,3];
    hits := [];
    print "N,cusp_form_dimension,full_cuspidal_symbols_dimension";
    for N in [1..LevelBound] do
        d := DimensionCuspFormsGamma1(N,3);
        symbols := Dimension(CuspidalSubspace(ModularSymbols(Gamma1(N),3,0)));
        assert symbols eq 2*d;
        if N le 12 then assert d eq small[N]; else assert d ge 5; end if;
        if d eq 2 then Append(~hits,N); end if;
        if N ge 5 then
            index_SL := N^2*&*[Rationals()|1-1/p^2 : p in PrimeDivisors(N)];
            cusps := &+[EulerPhi(a)*EulerPhi(N div a) : a in Divisors(N)]/2;
            // For N>=5 these groups are torsion-free with regular cusps.
            assert d eq index_SL/12-cusps/2;
        end if;
        printf "%o,%o,%o\n",N,d,symbols;
    end for;
    assert hits eq [9,10];
    expected_widths := AssociativeArray();
    expected_widths[9] := [1,1,1,3,3,9,9,9];
    expected_widths[10] := [1,1,2,2,5,5,10,10];
    expected_widths[12] := [1,1,2,3,3,4,4,6,12,12];
    for N in [9,10,12] do
        widths,mu := Gamma1CuspWidths(N);
        assert widths eq expected_widths[N];
        assert DimensionCuspFormsGamma1(N,3) eq (2*mu)/12-#widths/2;
        if N eq 12 then assert 2*mu eq 96 and #widths eq 10; end if;
        printf "Gamma1(%o): SL2 index=%o; cusps=%o; widths=%o; pg=%o\n",
               N,2*mu,#widths,widths,DimensionCuspFormsGamma1(N,3);
    end for;
    print "Quadratic refinements: N,discriminant,dimension,new_dimension";
    for pair in [<12,-4>,<14,-7>,<15,-3>,<15,-15>,<16,-8>,<18,-3>] do
        N,D := Explode(pair);
        chi := DirichletGroup(N,Rationals())!KroneckerCharacter(D);
        assert chi(-1) eq -1;
        S := CuspidalSubspace(ModularSymbols(chi,3,0));
        assert Dimension(S) eq 4;
        new_dim := Dimension(NewSubspace(S));
        assert new_dim eq (N in [14,16] select 0 else 4);
        printf "%o,%o,%o,%o\n",N,D,Dimension(S) div 2,new_dim div 2;
    end for;
    printf "PASS: Gamma1 dimensions through level %o; cusp widths for levels 9, 10, 12; six quadratic refinements.\n",LevelBound;
end procedure;

VerifyModularDimensions(40);
quit;
