// Exact Euler-polynomial utilities for rank-four weight-two motives.
// Convention: P_p(T) = det(1 - Frob_p*T), reciprocal roots of size p.

function OrthogonalEulerFromTraces(p, s1, s2, eta)
    assert IsPrime(p) and eta in {-1, 1};
    b := (s1^2-s2)/2;
    assert IsIntegral(b);
    assert eta eq 1 or b eq 0;
    R<T> := PolynomialRing(Integers());
    return 1-s1*T+(Integers()!b)*T^2-eta*p^2*s1*T^3+eta*p^4*T^4;
end function;

procedure CheckOrthogonalEuler(p, P, eta)
    assert IsPrime(p) and eta in {-1, 1};
    assert Degree(P) eq 4 and Coefficient(P,0) eq 1;
    assert Coefficient(P,4) eq eta*p^4;
    assert Coefficient(P,3) eq eta*p^2*Coefficient(P,1);
    assert eta eq 1 or Coefficient(P,2) eq 0;
end procedure;

function QMCandidateOrders(p, P)
    CheckOrthogonalEuler(p,P,1);
    a := -Coefficient(P,1); b := Coefficient(P,2);
    orders := {};
    if a^2 eq 4*b-8*p^2 then Include(~orders,1); end if;
    if a eq 0 then Include(~orders,2); end if;
    if a^2 eq b+p^2 then Include(~orders,3); end if;
    if a^2 eq 2*b then Include(~orders,4); end if;
    if a^2 eq 3*b-3*p^2 then Include(~orders,6); end if;
    return orders; // A set: special Frobenius polynomials can be ambiguous.
end function;

function TensorEulerPolynomial(a,d,b,e)
    // Factors have characteristic polynomials X^2-a*X+d and X^2-b*X+e.
    R<T> := PolynomialRing(Parent(a));
    return 1-a*b*T+(e*a^2+d*b^2-2*d*e)*T^2
             -a*b*d*e*T^3+d^2*e^2*T^4;
end function;

function AsaiInertEulerPolynomial(a,d)
    // Untwisted Asai+: a,d are trace and determinant at the inert prime of F.
    // Apply the weight-two normalization afterwards by T -> xi(p)*T/p.
    R<T> := PolynomialRing(Parent(a));
    return (1-a*T+d*T^2)*(1-d*T^2);
end function;
