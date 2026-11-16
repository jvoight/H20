// Rational quaternion conjugation matrices for the nine compact QM candidates.
// These are algebraic-group candidates, not Galois or surface realizations.
// Run from the 202-candidates directory: magma -b scripts/qm_candidates.m
load "scripts/euler_tools.m";

procedure VerifyQMCandidates()
    Q := Rationals();
    B<i,j,k> := QuaternionAlgebra<Q | -1,-1>;
    function Rotation(u)
        return Matrix(Q,3,3,&cat[Eltseq(u*v*u^-1)[2..4] : v in [i,j,k]]);
    end function;
    a3 := 1+i+j+k; a4 := 1+i; a6 := 3+i+j+k; s := i-j;
    names := ["C2","C4","C6","D2","D3","D4","D6","A4","S4"];
    generators := [ [i], [a4], [a6], [i,j], [a3,s], [a4,j],
                    [a6,s], [i,a3], [i,a3,a4] ];
    expected := [ [1,1,0,0,0], [1,1,0,2,0], [1,1,2,0,2],
                  [1,3,0,0,0], [1,3,2,0,0], [1,5,0,2,0],
                  [1,7,2,0,2], [1,3,8,0,0], [1,9,8,6,0] ];
    orders := [1,2,3,4,6];
    print "group,order,n1,n2,n3,n4,n6,trace_zero_density";
    for r in [1..#names] do
        mats := [Rotation(u) : u in generators[r]];
        for A in mats do
            assert Determinant(A) eq 1 and A*Transpose(A) eq IdentityMatrix(Q,3);
        end for;
        G := sub<GL(3,Q) | mats>;
        counts := [#[g : g in G | Order(g) eq n] : n in orders];
        assert counts eq expected[r] and &+counts eq #G;
        printf "%o,%o,%o,%o,%o,%o,%o,%o\n", names[r],#G,
                counts[1],counts[2],counts[3],counts[4],counts[5],Q!counts[2]/#G;
    end for;
    // Verify the tensor-polynomial formula symbolically, not just numerically.
    R<a,b,d,e> := PolynomialRing(Q,4);
    A := Matrix(R,2,2,[0,-d,1,a]); C := Matrix(R,2,2,[0,-e,1,b]);
    P := CharacteristicPolynomial(TensorProduct(A,C));
    assert [Coefficient(P,4-r):r in [0..4]] eq
       [R|1,-a*b,e*a^2+d*b^2-2*d*e,-a*b*d*e,d^2*e^2];
    assert [Coefficient(P,4-r):r in [0..4]] eq
        Coefficients(TensorEulerPolynomial(a,d,b,e));
    Swap := Matrix(R,4,4,[1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1]);
    PA := CharacteristicPolynomial(TensorProduct(A,IdentityMatrix(R,2))*Swap);
    assert [Coefficient(PA,4-r):r in [0..4]] eq
        Coefficients(AsaiInertEulerPolynomial(a,d));
    P0 := OrthogonalEulerFromTraces(5,6,-38,1);
    assert Coefficients(P0) eq [1,-6,37,-150,625];
    assert 6 in QMCandidateOrders(5,P0);
    P1 := OrthogonalEulerFromTraces(5,3,9,-1);
    CheckOrthogonalEuler(5,P1,-1);
    print "PASS: nine rational rotation groups; tensor, Asai, and trace identities.";
end procedure;

VerifyQMCandidates();
quit;
