MaxMultiplicity := function(S);
  m := 1;
  for x in S do
    m0 := Multiplicity(S,x);
    if m0 gt m then
      t := x;
      m := m0;
    end if;
  end for;    
  return t, m;
end function;  

for ell := 1 to #lfss do
  print "===============", ell;
  Fs := {* *};
  for i := 1 to #lfss[ell] do
	f := lfss[ell][i][2];
	K := SplittingField(f);
	if K cmpeq Rationals() then continue; end if;
	gs := {* Polredabs(DefiningPolynomial(Fdat[1])) : Fdat in Subfields(K) *};
	Fs join:= gs;
  end for;
  f, m := MaxMultiplicity(Fs);
  f, m;
  F := NumberField(f);
  [Factorization(f[2],F) : f in lfss[ell]];  
end for;

BaseChange := function(f,r);
  M := MatrixWithGivenCharacteristicPolynomial(Reverse(f));
  return Reverse(CharacteristicPolynomial(M^r));
end function;

for lfd in lfs do
  p := lfd[1];
  lf := lfd[2];
  if LegendreSymbol(-5,p) eq -1 then
    lf := BaseChange(lf,2);
  end if;
  print p, LegendreSymbol(-5,p), Factorization(lf,CyclotomicField(4));
end for;

psi := HeckeCharacterGroup(2*Integers(K)).1;
Kz := Codomain(psi);
_<xz> := PolynomialRing(Kz);
for lfd in lfs do
  p := lfd[1];
  lf := lfd[2];
  lf := Evaluate(lf,xz/p);
  if Degree(lf) lt 4 then continue; end if;
  lffact := Factorization(lf);
  if LegendreSymbol(-5,p) ne 1 then continue; end if;
  if #lffact eq 2 and &and[Degree(lfi[1]) eq 2 : lfi in lffact] then
    pps := Factorization(p*Integers(K));
    assert {* Coefficient(lfFacti[1],0)^^lfFacti[2] : lfFacti in lffact *} eq 
           {* psi(ppi[1]) : ppi in pps *};
    if Evaluate(lffact[1],xz*psi(pps[1][1])) 
  end if;
  print p, LegendreSymbol(-5,p), Factorization(lf,Kz), {* psi(ppi[1]) : ppi in Factorization(p*Integers(K)) *};
end for;

CM or partly CM 
17 : QQ(i)    640.3.e.c?
30 : QQ(sqrt(8),sqrt(-3))
35 : QQ(i)  

reducible over QQ
1 : CMF, coeff QQ(i)    32.3.c.a
2 : CMF, coeff QQ(omega)
3 : CMF, coeff QQ(sqrt(-8))   24.3.e.a
4 : CMF, coeff QQ(i)
13, 15 : CMF, coeff QQ(zeta_8)     13 and 15 quadratic twists (3-isogeny)
16 : CMF, coeff QQ(sqrt(-24))
23 : CMF, coeff QQ(i)    32.3.c.a
24 : CMF, coeff QQ(zeta_8)
25 : CMF, coeff QQ(i)
28 : CMF, coeff QQ(i)    288.3.g.a
29 : CMF, coeff QQ(zeta_8)
32 : CMF, coeff QQ(zeta_12,sqrt(-24))

potential reducible over QQ
14 : CMF, coeff QQ(zeta_8)     13 and 15 quadratic twists (3-isogeny)
18 : CMF, coeff QQ(zeta_12)

irreducible
5-8 : Asai HMF over QQ(sqrt(12))
9-12 : Asai HMF over QQ(sqrt(8))    10 and 12 quadratic twists (2-isogeny)
19-20 : Asai HMF over QQ(sqrt(24))   
21-22 : Asai HMF over QQ(sqrt(12))
26-27 : Asai HMF over QQ(sqrt(8))
31 : Asai HMF over QQ(sqrt(24))
33 : Asai HMF over QQ(sqrt(12))
34 : Asai HMF over QQ(sqrt(8))

