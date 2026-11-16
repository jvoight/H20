// Run from verification/: magma -b herfurtner_checks.m
// Reconstruction of the checks cited in rewrite-v11.tex, not the lost source.
// Recompute geometry and complete L-polynomials; manuscript values are oracles.
load "herfurtner_models.m";
load "../202-candidates/scripts/euler_tools.m";

function GeometricFibers(E)
    return Sort(&cat[[Sprint(k[1]) : j in [1..k[2]]] : k in KodairaSymbols(E)]);
end function;

function RootDiscriminant(k)
    if k eq KodairaSymbol("II*") then return 1; end if;
    if k eq KodairaSymbol("III*") then return 2; end if;
    if k eq KodairaSymbol("IV*") then return 3; end if;
    assert k eq KodairaSymbol("In*");
    return 4;
end function;

procedure VerifyHerfurtner()
    SetSeed(1);
    SetColumns(0);
    v1,v2,v3 := GetVersion();
    printf "Herfurtner checks, Magma %o.%o-%o\n",v1,v2,v3;
    primes := [5,7,13,17,19,31,37];
    // Rows of v11's full-Euler diagnostics: traces, quadratic coefficients,
    // and the complete sets of compatible projective orders (singletons here).
    expected := AssociativeArray();
    expected[13] := <[-6,0,0,18,0,0,68],[18,-82,-274,162,-466,14,3894],
                     [4,2,2,4,2,2,1]>;
    expected[15] := <[6,0,0,-18,0,0,68],[18,-82,-274,162,-466,14,3894],
                     [4,2,2,4,2,2,1]>;
    expected[25] := <[0,0,4,0,0,0,36],[-34,-82,342,-2,-146,782,3062],
                     [2,2,1,2,2,2,1]>;
    expected[28] := <[8,0,-4,48,0,0,36],[66,-82,342,1154,-146,782,3062],
                     [1,2,1,1,2,2,1]>;
    h7 := AssociativeArray();
    h7[13] := [1,-4,-186,-676,28561];
    h7[37] := [1,-52,1158,-71188,1874161];
    fibers := AssociativeArray();
    fibers[7] := ["IV*","I5*","I4*","I1*"];
    fibers[13] := ["II*","I6*","I1*","I1*"];
    fibers[15] := ["II*","I3*","I3*","I2*"];
    fibers[25] := ["III*","III*","I3*","I3*"];
    fibers[28] := ["II*","II*","I2*","I2*"];
    ZT<T> := PolynomialRing(Integers());
    total := 0; diagnostic_rows := 0; d4_checks := 0;
    print "Coefficient rows: model,p,c0,c1,c2,c3,c4,eta";
    for id in [7,13,15,25,28] do
        E,d := HerfurtnerInput(id);
        Q := BaseRing(E);
        // Check the explicit twist against the original construction.
        assert d eq &*[Numerator(Q!Minimum(k[1])) : k in LocalInformation(E)];
        X := QuadraticTwist(E,d);
        assert GeometricFibers(X) eq Sort(fibers[id]);
        assert BettiNumber(X,1) eq 0 and BettiNumber(X,2) eq 34;
        ks := KodairaSymbols(X);
        trivial_rank := 2+&+[k[2]*(NumberOfComponents(k[1])-1) : k in ks];
        // e=36, chi(O)=3, pg=2, h11=30. Shioda--Tate gives MW rank zero.
        assert trivial_rank eq 30;
        root_disc := &*[RootDiscriminant(k[1])^k[2] : k in ks];
        orientation_disc := id eq 7 select 3 else 1;
        assert IsSquare(Rationals()!root_disc/orientation_disc);
        printf "H%o: twist=%o; geometric fibers=%o; b2=34; trivial rank=30; orientation squareclass=%o\n",
               id,d,GeometricFibers(X),orientation_disc;
        // H7 uses the two primes displayed in v11. In particular, p=5 is
        // bad for this model: its finite bad places t=0 and t=-10 collide.
        model_primes := id eq 7 select [13,37] else primes;
        for j in [1..#model_primes] do
            p := model_primes[j];
            Xp := ChangeRing(X,FunctionField(GF(p)));
            assert GeometricFibers(Xp) eq GeometricFibers(X);
            assert BettiNumber(Xp,2) eq 34;
            // LFunction removes the trivial lattice, which is all of the
            // characteristic-zero NS contribution. Retain every other factor,
            // including factors from new algebraic classes after reduction.
            P := ZT!LFunction(Xp);
            eta := KroneckerSymbol(orientation_disc,p);
            CheckOrthogonalEuler(p,P,eta);
            C := Coefficients(P);
            printf "H%o,%o,%o,%o,%o,%o,%o,%o\n",
                   id,p,C[1],C[2],C[3],C[4],C[5],eta;
            if id eq 7 then
                if eta eq -1 then
                    assert P eq (1-p^2*T^2)*(1+C[2]*T+p^2*T^2);
                end if;
                if IsDefined(h7,p) then
                    assert C eq h7[p];
                    assert IsIrreducible(P);
                    G := GaloisGroup(P);
                    // An irreducible quartic with Galois group of order 8
                    // has the dihedral group D4 (of order 8).
                    assert #G eq 8 and not IsAbelian(G);
                    printf "H7,p=%o: irreducible quartic; Galois group D4 (order %o)\n",p,#G;
                    d4_checks +:= 1;
                end if;
            else
                a := expected[id][1][j]; b := expected[id][2][j];
                assert C eq [1,-a,b,-p^2*a,p^4];
                orders := QMCandidateOrders(p,P);
                assert orders eq {expected[id][3][j]};
                printf "H%o,p=%o: (a,b)=(%o,%o); compatible orders=%o\n",id,p,a,b,orders;
                diagnostic_rows +:= 1;
            end if;
            total +:= 1;
        end for;
    end for;
    assert total eq 30 and diagnostic_rows eq 28 and d4_checks eq 2;
    print "PASS: five characteristic-zero models; 30 full Euler polynomials; 28 manuscript diagnostic entries; two H7 D4 quartics.";
end procedure;

VerifyHerfurtner();
quit;
