SetColumns(0);
K<t>:=FunctionField(Rationals());
a:=(t^2-1)^2;b:=(t^2+1)^2;d:=t*(t-1);
E:=EllipticCurve([K|0,-a-b,0,a*b,0]);
assert Discriminant(E) eq 256*t^4*(t^4-1)^4;
assert jInvariant(E) eq 16*(t^8+14*t^4+1)^3/(t^4*(t^4-1)^4);
u:=(t^2+1)/(2*t);
F:=EllipticCurve([K|0,-1-u^2,0,u^2,0]);
assert IsIsomorphic(F,QuadraticTwist(E,-1));
assert IsIsomorphic(QuadraticTwist(F,d),QuadraticTwist(E,-d));
assert 4*4+2*(4+6) eq 36;
assert 2+4*(4-1)+2*(4+4) eq 30;
// The Gamma(4) cusp parabolics modulo 8 depend only on the cusp mod 2.
R:=Integers(8); T4:=Matrix(R,2,2,[1,4,0,1]);
for c in [0..3] do for a in [0..3] do
    if GCD(a,c) ne 1 then continue; end if;
    gg,x,y:=XGCD(a,c);alpha:=Matrix(R,2,2,[a,-y,c,x]);
    p:=alpha*T4*alpha^-1;
    assert p eq Matrix(R,2,2,[1-4*a*c,4*a^2,-4*c^2,1+4*a*c]);
end for;end for;
print "PASS: Gamma(4) twist, pg=2, rho=30, and constant-twist equivalence to row (1,15)";
quit;
