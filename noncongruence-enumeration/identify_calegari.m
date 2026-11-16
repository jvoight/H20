// Numerical monodromy, independently repeated with a finer path below.
SetColumns(0);
CC<I>:=ComplexField(80); RR:=RealField(80); P<u>:=PolynomialRing(CC);
phi:=3*u^2-2*u^3;
C:=(phi+2)*(phi^3+6*phi^2-12*phi+8);
F:=C^3; D:=1728*phi^3*(phi-1)^2*(phi+8);
Fp:=Derivative(F); Dp:=Derivative(D);
roots:=[r[1]:r in Roots(F-D/2)];assert #roots eq 36;
function Monodromy(center,steps)
    values:=roots; previous:=(CC!1)/2;
    for k in [1..steps] do
        z:=center+((CC!1)/2-center)*Exp(2*Pi(RR)*I*k/steps);
        next:=[];
        for x in values do
            y:=x+(z-previous)*Evaluate(D,x)/(Evaluate(Fp,x)-previous*Evaluate(Dp,x));
            for it in [1..20] do
                delta:=(Evaluate(F,y)-z*Evaluate(D,y))/(Evaluate(Fp,y)-z*Evaluate(Dp,y));
                y-:=delta;if Abs(delta) lt 10^-55 then break;end if;
            end for;
            assert Abs(Evaluate(F,y)-z*Evaluate(D,y))/(1+Abs(Evaluate(F,y))) lt 10^-50;
            Append(~next,y);
        end for;
        assert Min([Abs(next[i]-next[j]):i,j in [1..36]|i lt j]) gt 10^-20;
        values:=next; previous:=z;
    end for;
    perm:=[];
    for y in values do
        ds:=[Abs(y-x):x in roots]; v,j:=Min(ds);assert v lt 10^-40;Append(~perm,j);
    end for;
    return Sym(36)!perm;
end function;
r0:=Monodromy(0,240);s:=Monodromy(1,240);
assert r0 eq Monodromy(0,480) and s eq Monodromy(1,480);
r:=r0^-1;
assert Order(r) eq 3 and Order(s) eq 2;
assert IsTransitive(sub<Sym(36)|s,r>);
assert Sort([#o:o in Orbits(sub<Sym(36)|s*r>)]) eq [1,1,1,2,3,4,6,18];
PrintFile("calegari-action.txt",Sprintf("%o\n%o",[i^r-1:i in [1..36]],[i^s-1:i in [1..36]]):Overwrite:=true);
print "PASS: repeated numerical monodromy and exact permutation passport";
quit;
