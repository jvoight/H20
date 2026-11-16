// Run from the 202-candidates directory: magma -b scripts/example_S1_9.m
// All primes p <= 101 except p=3. No online data or precomputed eigenvalues.
// For a different bound edit PrimeBound below; the q-expansion grows with it.
load "scripts/euler_tools.m";

procedure RunS19(PrimeBound)
    v1,v2,v3 := GetVersion();
    printf "S1(9), Magma version: %o.%o-%o; bound: %o\n", v1,v2,v3,PrimeBound;
    assert PrimeBound ge 2;
    K<z> := CyclotomicField(6);
    D := DirichletGroup(9,K); chi := D.1;
    assert Order(chi) eq 6 and chi(-1) eq -1;
    S := CuspidalSubspace(ModularSymbols(chi,3,1));
    assert Dimension(S) eq 1;
    f := qEigenform(S,PrimeBound+1);
    sigma := hom<K->K | 1-z>;
    ZT<T> := PolynomialRing(Integers());
    KT<U> := PolynomialRing(K);

    // The characteristic-zero elliptic surface has fibers 3 I9, 2 I3, 3 I1.
    // Euler number 36 gives chi(O)=3, pg=2, b2=34, h11=30.
    // Its trivial lattice has rank 2+3*8+2*2=30, hence geometric MW rank 0.
    // Therefore LFunction(E_p) removes exactly the specialization of NS(X_Qbar)
    // at the good primes used here. Do NOT remove additional factors coming
    // from Picard-rank jumps in characteristic p.
    Q<t> := FunctionField(Rationals());
    c := t^2*(t-1); b := c*(t^2-t+1);
    assert b^3-3*b^2*c+b*c^3+3*b*c^2-c^5-c^4-c^3 eq 0;
    E := EllipticCurve([Q | 1-c,-b,-b,0,0]);
    printf "Characteristic-zero fibers: %o\n", KodairaSymbols(E);
    printf "chi(2) = %o; f = %o\n",chi(2),qEigenform(S,12);
    print "p,projective_order,c0,c1,c2,c3,c4";
    checked := [];
    for p in PrimesUpTo(PrimeBound) do
        if p eq 3 then continue; end if;
        F<t> := FunctionField(GF(p));
        c := t^2*(t-1); b := c*(t^2-t+1);
        Ep := EllipticCurve([F | 1-c,-b,-b,0,0]);
        Pgeom := ZT!LFunction(Ep);
        a := K!Coefficient(f,p); d := chi(p)*p^2;
        Pmod := ZT!((1-a*U+d*U^2)*(1-sigma(a)*U+sigma(d)*U^2));
        assert sigma(a) eq a/chi(p);
        // Artin factor tau = 1 direct-sum chi^(-1), no additional twist.
        Ptensor := ZT!TensorEulerPolynomial(a,d,1+1/chi(p),1/chi(p));
        assert Pgeom eq Pmod and Pmod eq Ptensor;
        CheckOrthogonalEuler(p,Pgeom,1);
        n := Min([n : n in [1,2,3,6] | chi(p)^n eq 1]);
        assert n in QMCandidateOrders(p,Pgeom);
        C := Coefficients(Pgeom);
        printf "%o,%o,%o,%o,%o,%o,%o\n",p,n,C[1],C[2],C[3],C[4],C[5];
        Append(~checked,p);
    end for;
    assert checked eq [p : p in PrimesUpTo(PrimeBound) | p ne 3];
    printf "PASS: %o full Euler polynomials; primes <= %o, excluding 3.\n",
             #checked,PrimeBound;
end procedure;

RunS19(101);
quit;
