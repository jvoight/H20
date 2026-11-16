Hs := [HypergeometricData(H) : H in PossibleHypergeometricData(5 : Weight := 2) |
       HodgeVector(HodgeStructure(HypergeometricData(H),1)) eq [2,0,2]];
_<T> := PolynomialRing(Integers());

K := QuadraticField(-1);

for H in [Hs[i] : i in [3, 7, 13, 14, 19, 27]] do
  LH1 := LSeries(H,1);
  print "========", Index(Hs, H);
  print H;
  print [<p,Factorization(EulerFactor(LH1,p),K)> : p in PrimesUpTo(200)];
  
  [[Squarefree(Discriminant(K[1])) : K in Subfields(SplittingField(EulerFactor(LH1,p))) | 
    Degree(K[1]) eq 2] : p in PrimesUpTo(200) | Degree(EulerFactor(LH1,p)) gt 1 
    and Degree(SplittingField(EulerFactor(LH1,p))) gt 1];
end for;

[ 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 ]

p := 11;
_<T> := PolynomialRing(Rationals());
// Lp := (1-30*T+p^2*T^2)*(1+16*T+p^2*T^2);
// Lp := 1-335*T^2+p^4*T^4;
Lp := (1-10*T+p^2*T^2)^2;
L17s := [];
for i in [1..#Hs] do
  H := Hs[i];
  for t in [1..p-1] do
    LHtp := EulerFactor(LSeries(H,t),p);
    if Degree(LHtp) eq 5 then
      if LHtp mod (1-p*T) eq 0 then 
        Append(~L17s, <i,t,LHtp div (1-p*T)>);
      end if;
      if LHtp mod (1+p*T) eq 0 then 
        Append(~L17s, <i,t,LHtp div (1+p*T)>);
      end if;
    else
      if Degree(LHtp) eq 4 then
        Append(~L17s, <i,t,LHtp>);
      end if;
    end if;
  end for;
end for;

allpds := [];
for H in Hs do
  pds := AssociativeArray();
  bads := [];
  for p in PrimesInInterval(11,50) do
	print p;
	for t in [1..p-1] do
	  LHt := LSeries(H,t);
	  LHt_spl := SplittingField(EulerFactor(LHt,p));
	  try
		if LHt_spl eq Rationals() then
		  ds := [1];
		else  
		  ds := [Squarefree(Discriminant(K[1])) : 
					 K in Subfields(LHt_spl) | 
					 Degree(K[1]) le 2];
		end if;
		print t, ds;
		for d in ds do
		  if d in Keys(pds) then
			Append(~pds[d], <p,t>);
		  else
			pds[d] := [<p,t>];
		  end if;
		end for;
	  catch e;
		Append(bads, <p,t,EulerFactor(LHt,p)>);
	  end try;
	end for;
  end for;
  for d in Keys(pds) do
	if d gt 0 or Abs(d) gt 50 then continue; end if;
	print d, pds[d];
  end for;
  Append(~allpds, pds);
end for;
