// 289.3.e.o in a rational rank-four QM tensor product.
// Run from verification/: MAGMA_ARCH=arm64 magma -b qm_289_counterexample.m
// No downloaded coefficients or external packages are used.
// See README.md for references, the mathematical descent argument, and scope.
// All checks lie in one procedure: an assertion failure prevents the final PASS.

function RamificationSet(B)
    S := Seqset(PrimeDivisors(Integers()!Discriminant(B)));
    if IsDefinite(B) then Include(~S,0); end if; // 0 denotes infinity.
    return S;
end function;

function AddRamification(S,T)
    return (S diff T) join (T diff S);
end function;

procedure VerifyQM289(PrimeBound)
    assert PrimeBound ge 89;
    SetSeed(1);
    SetColumns(0);
    v1,v2,v3 := GetVersion();
    printf "QM 289.3.e.o verification, Magma %o.%o-%o\n",v1,v2,v3;

    // Unique degree-16 orbit in S_3(289,chi), chi primitive of conductor 17.
    G := FullDirichletGroup(289);
    chars := [c : c in GaloisConjugacyRepresentatives(G) |
                 IsOdd(c) and Order(c) eq 16];
    assert #chars eq 1;
    chi := MinimalBaseRingCharacter(chars[1]);
    C := BaseRing(chi); z := C.1;
    powers := [j : j in [1..15 by 2] | Evaluate(chi^j,3) eq z];
    assert #powers eq 1;
    chi := chi^powers[1];
    assert Conductor(chi) eq 17 and Order(chi) eq 16;
    assert Evaluate(chi,3) eq z and Evaluate(chi,2) eq z^14;
    pieces := NewformDecomposition(NewSubspace(
                  CuspidalSubspace(ModularSymbols(chi,3,1))));
    selected := [A : A in pieces | Dimension(A) eq 2];
    assert #selected eq 1;
    prec := Maximum(90,PrimeBound+1);
    f := qEigenform(selected[1],prec);
    K := BaseRing(Parent(f));
    assert AbsoluteDegree(K) eq 16;

    // Fix a presentation independent of the generator returned by qEigenform.
    PC<X> := PolynomialRing(C);
    E<b> := NumberField(X^2-22*z);
    ze := E!z;
    roots := Roots(DefiningPolynomial(K),E);
    assert #roots eq 2;
    embeddings := [hom<K -> E | root[1]> : root in roots];
    embeddings := [h : h in embeddings | h(Coefficient(f,3)) eq b];
    assert #embeddings eq 1;
    hf := embeddings[1];
    coeff := [hf(Coefficient(f,n)) : n in [0..prec-1]];
    assert coeff[3] eq ze^7 and coeff[4] eq b;
    assert coeff[3]^7 eq ze; // a_2 and a_3 generate E, not just a subfield.
    assert &and[coeff[n+1] eq 0 : n in [1..prec-1] | n mod 17 eq 0];
    print "Field: E = Q(zeta_16,b), b^2 = 22*zeta_16; [E:Q] = 16.";
    print "Nebentypus: conductor 17, chi(3)=zeta_16; a_2=zeta_16^7, a_3=b.";

    // Twisting by chi^j has level dividing lcm(289,17^2,17*17)=289.
    // The common-character Sturm bound is floor(3*[SL2:Gamma0(289)]/12)=76.
    sturm := Floor(3*289*(1+Rationals()!1/17)/12);
    assert sturm eq 76;
    print "Exact inner twists: r, sign on b, exponent j in f^sigma=f tensor chi^j";
    for r in [1..15 by 2] do
        for s in [-1,1] do
            hbase := hom<C -> E | ze^r>;
            sigma := hom<E -> E | hbase,s*ze^((r-1) div 2)*b>;
            j := (r-1) div 2 + (s eq -1 select 8 else 0);
            assert &and[sigma(coeff[n+1]) eq
                         (E!Evaluate(chi^j,n))*coeff[n+1] : n in [1..89]];
            assert r mod 16 eq (1+2*j) mod 16;
            printf "%o, %o, %o\n",r,s,j;
        end for;
    end for;
    print "PASS: all 16 inner twists through q^89 (Sturm bound 76).";
    // No imaginary quadratic fundamental discriminant divides 289, so f is non-CM.

    // The scalar alpha(m) and its rational multiplication cocycle.
    alpha := [m mod 2 eq 0 select ze^(m div 2)
              else b*ze^((m-1) div 2) : m in [0..15]];
    for m,n in [0..15] do
        c := alpha[m+1]*alpha[n+1]/alpha[((m+n) mod 16)+1];
        expected := 22^((m mod 2)*(n mod 2))*(-1)^((m+n) div 16);
        assert c eq expected;
    end for;

    // Brauer classes: gamma_chi has ramification {17,infinity}, since the
    // local character is odd at 17 and chi is odd globally. The positive
    // part of c_alpha contributes (17,22). No numerical InnerTwists call.
    gamma := {0,17};
    B<I,J,IJ> := QuaternionAlgebra<Rationals() | -22,-1>;
    assert I^2 eq -22 and J^2 eq -1 and I*J eq -J*I;
    positive := QuaternionAlgebra<Rationals() | 17,22>;
    assert AddRamification(gamma,RamificationSet(positive)) eq RamificationSet(B);
    assert RamificationSet(B) eq {0,11};
    sign_symbols := [[17,-1],[-3,-1],[17,-3]];
    sign_ramification := {};
    for ab in sign_symbols do
        A := QuaternionAlgebra<Rationals() | ab[1],ab[2]>;
        sign_ramification := AddRamification(sign_ramification,RamificationSet(A));
        printf "Sign symbol %o: ramification %o (0=infinity)\n",ab,Sort(Setseq(RamificationSet(A)));
    end for;
    assert sign_ramification eq gamma;
    for e,h,ee,hh in [0..1] do
        q1 := I^e*J^h; q2 := I^ee*J^hh;
        c := 22^(e*ee)*(-1)^(e*ee+h*hh+h*ee);
        assert q1*q2 eq c*I^((e+ee) mod 2)*J^((h+hh) mod 2);
        assert Norm(q1) eq 22^e;
    end for;
    print "PASS: B has discriminant 11; the two sign cocycles have equal Brauer class.";

    // Concrete odd Artin factor: induce an order-32 ray character from Q(sqrt(-3)).
    L<w> := QuadraticField(-3); O := MaximalOrder(L);
    P3 := Factorisation(3*O)[1][1];
    R,mp := RayClassGroup(17*O*P3);
    R17 := RayClassGroup(17*O);
    R3 := RayClassGroup(P3);
    assert ClassNumber(L) eq 1 and Invariants(R) eq [96];
    assert #R17 eq 48 and #R3 eq 1;
    // An order-32 character cannot factor through either proper divisor:
    // its conductor is exactly 17*P3. Induction has conductor 3*Norm(17*P3)=2601.
    C32<t> := CyclotomicField(32);
    embedC := hom<C -> C32 | t^2>;
    units51 := [n : n in [1..50] | GCD(n,51) eq 1];
    good_powers := [u : u in [1..31 by 2] |
        &and[KroneckerSymbol(-3,n)*t^(u*Eltseq((n*O)@@mp)[1]) eq
             embedC(Evaluate(chi,n))^-1 : n in units51]];
    assert #good_powers eq 2;
    u := good_powers[1];
    assert good_powers[2] eq (17*u) mod 32;
    conj := hom<L -> L | -w>;
    ray_generator := mp(R.1);
    conjugate_generator := ideal<O | [conj(L!a) : a in Basis(ray_generator)]>;
    nrm := Integers()!Norm(ray_generator);
    assert GCD(nrm,51) eq 1;
    ratio := t^(u*(Eltseq(conjugate_generator@@mp)[1]-Eltseq(ray_generator@@mp)[1]));
    assert ratio eq KroneckerSymbol(17,nrm);
    // This generator check proves theta^conj/theta=chi_17 restricted to G_L.
    printf "Ray group C_96; compatible order-32 character powers %o.\n",good_powers;
    print "PASS: Artin determinant chi^-1; conductor 2601; projective field Q(sqrt(17),sqrt(-3)).";

    // Exact Euler-polynomial cross-checks; these are finite checks, not the
    // proof of descent. The latter uses B_l tensor B_l = M_4(Q_l).
    PT<Y> := PolynomialRing(C32); M<s> := NumberField(Y^2-22);
    mapC := hom<C -> M | (M!t)^2>;
    mapE := hom<E -> M | mapC,s*(M!t)>;
    PM<T> := PolynomialRing(M);
    primes := [p : p in PrimesInInterval(2,PrimeBound) | p notin [3,17]];
    print "p, trace, second coefficient, Euler coefficients (ascending powers)";
    for p in primes do
        pp := Factorisation(p*O);
        trtau := #pp eq 1 select C32!0
                 else &+[t^(u*Eltseq(P[1]@@mp)[1]) : P in pp];
        dettau := embedC(Evaluate(chi,p))^-1;
        r := trtau^2/dettau;
        assert r in {C32!0,C32!4};
        assert (r eq 4) eq
               (KroneckerSymbol(17,p) eq 1 and KroneckerSymbol(-3,p) eq 1);
        ap := mapE(coeff[p+1]);
        chip := M!embedC(Evaluate(chi,p));
        mu := ap^2/chip;
        m := [j : j in [0..15] | ze^j eq E!Evaluate(chi,p)][1];
        assert IsCoercible(Rationals(),coeff[p+1]/alpha[m+1]);
        assert IsCoercible(Rationals(),mu);
        aT := ap*(M!trtau);
        bT := mu+(M!r-2)*p^2;
        assert IsCoercible(Integers(),aT) and IsCoercible(Integers(),bT);
        expected := [M!1,-aT,bT,-p^2*aT,M!p^4];
        rfrob := Matrix(M,2,2,[0,-p^2*chip,1,ap]);
        tfrob := Matrix(M,2,2,[0,-M!dettau,1,M!trtau]);
        cp := CharacteristicPolynomial(TensorProduct(rfrob,tfrob));
        actual := [Coefficient(cp,4-j) : j in [0..4]];
        assert actual eq expected;
        printf "%o, %o, %o, %o\n",p,Integers()!aT,Integers()!bT,
               [Integers()!x : x in actual];
    end for;
    printf "PASS: %o complete rational integral tensor Euler polynomials at good primes <= %o.\n",#primes,PrimeBound;

    // Restrict an arbitrary Dirichlet twist to mu_16 in Z_17^times.
    // Squaring adds an even exponent, so the order-16 obstruction persists.
    assert &and[(z^(m+2*j))^16 eq 1 and (z^(m+2*j))^8 ne 1
                : m in [1..15 by 2], j in [0..15]];
    print "PASS: all possible local square twists preserve order 16 at 17.";
    print "PASS: qm_289_counterexample (exact inner twists, Brauer symbols, ray character, Euler factors).";
end procedure;

VerifyQM289(199);
quit;
