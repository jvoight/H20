P3w<x,y,z,w> := ProjectiveSpace(Q,[1,1,1,4]);

P2<a,b,c> := ProjectiveSpace(Q,2);
to2 := map<P3w->P2|[x,y,z]>;
// cf. Persson
steiner := (3*b^2+a^2)^2+b*c*(12*a^2-28*b^2-12*c^2)+c^2*(-2*a^2+30*b^2+c^2);
cs := Curve(P2,steiner);
cusps := SingularPoints(cs);
assert forall{c: c in cusps|IsCusp(c)};
bit := Curve(P2,c); // bitangent
tans := [Curve(PrimeComponents(TangentCone(i))[1]): i in cusps]; // tangents

// Let's find the automorphisms of the Steiner quartic.
pl := Places(cusps[1])[1];
dm := DivisorMap(1*pl);
_,md := IsInvertible(dm);
P1 := Codomain(dm);
ptsp1 := [dm(i): i in cusps];
trs := [TranslationOfSimplex(P1,[ptsp1[i^p]: i in [1..3]]): p in SymmetricGroup(3)|Order(p) ne 3];
tr2 := [trs[1]^-1*i: i in trs];
auts := [dm*i*md: i in tr2];
aut_small := [map<cs->cs|AllDefiningPolynomials(Extend(Normalization(Expand(i))))[2]>: i in auts];
assert forall{x: x in aut_small|Degree(DefiningPolynomials(x)[1]) eq 1};

// so now example 3.2.1 is the Steiner quartic together with these four lines.
ex321 := Surface(P3w,w^2-DefiningPolynomial(&join (tans cat [bit,cs]))@@to2);
// This has an action of a group of order 12 over Q(i): we can permute the
// three cusp tangents and change the sign of w.

// But if we apply an odd permutation, that takes w to iw, so the automorphism
// is not defined over Q.
Qi<i> := QuadraticField(-1);
ex321i := ChangeRing(ex321,Qi);
aut_321i := [map<ex321i->ex321i|Append([CoordinateRing(AmbientSpace(ex321i))|i@@to2: i in DefiningPolynomials(a)],i*ex321i.4)>: a in aut_small[2..4]];

aut_321 := [map<ex321->ex321|DefiningPolynomials(&*[aut_321i[i]: i in Remove([1..3],j)])>: j in [1..3]];

// This has a K3 surface quotient, if we mod out by an automorphism of order 3.
// The quotient has to do with
// http://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/36/3/d/a/.

// There is no other subgroup of this group of order 12 that gives a quotient
// with h^{2,0} > 0, so we don't necessarily get whatever it is on a surface.
// This seems to realize the representation attached to
// http://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/36/3/d/c/.
// counts are
// [ <5, 12>, <7, 0>, <11, 0>, <13, -6>, <17, -36>, <19, 0>, <23, 0>, <29, 12>, <31, 0>, <37, -18>, <41, -36>, <43, 0>, <47, 0>, <53, 204>, <59, 0>, <61, 30>, <67, 0>, <71, 0>, <73, 18>, <79, 0>, <83, 0>, <89, -324>, <97, -126> ]
// but we need to subtract off the quotient to see this.

// example 3.2.2
ls2 := LinearSystem(P2,2);
ls2 := LinearSystem(ls2,ReducedSubscheme(cs meet bit));
for i in [1,3] do ls2 := LinearSystem(ls2,P2!cusps[i]); end for;
// we can find it using infinitely near points or otherwise
con := Curve(P2,ls2.1+6*ls2.2);
ex322 := Surface(P3w,w^2-DefiningPolynomial(&join ([bit,tans[2],con,cs]))@@to2);
br322 := DefiningPolynomial(to2(Scheme(ex322,w)));
cts322 := [<p,p le 3 select 1 else &+[KroneckerSymbol(Evaluate(br322,Eltseq(i))): i in RationalPoints(ProjectiveSpace(GF(p),2))]>: p in PrimesUpTo(100)];

// this has an automorphism (x:y:z:w) -> (-x:y:z:iw)
// but the quotient will go through P2, so we don't decompose.

// example 3.2.3
// find using infinitely near points again
// this is for the cusp at (0:1:1)
con := Curve(P2,a^2+12*b^2+8*b*c-28/27*c^2);
ex323 := Surface(P3w,w^2-DefiningPolynomial(&join ([bit,tans[2],con,cs]))@@to2);
// again this has an automorphism (x:y:z:w) -> (-x:y:z:iw).
br323 := DefiningPolynomial(to2(Scheme(ex323,w)));
cts323 := [<p,p le 3 select 1 else &+[KroneckerSymbol(Evaluate(br323,Eltseq(i))): i in RationalPoints(ProjectiveSpace(GF(p),2))]>: p in PrimesUpTo(100)];

ex324 := Surface(P3w,w^2-DefiningPolynomial(&join ([tans[1],tans[2],con,cs]))@@to2);
br324 := DefiningPolynomial(to2(Scheme(ex324,w)));
cts324 := [<p,p le 3 select 1 else &+[KroneckerSymbol(Evaluate(br324,Eltseq(i))): i in RationalPoints(ProjectiveSpace(GF(p),2))]>: p in PrimesUpTo(100)];

// now let's look at the ad hoc examples
// example 3.3.1
// choosing the nodal cubic to make the nodal tangent rational
load "k3pic.mag";
nod := Curve(P2,b^2*c-a^2*(a+c));
fl := Flexes(nod);                
fl := excise(fl,Scheme(P2,[a,b]));
Qrm3<rm3> := QuadraticField(-3);
tans3 := [TangentLine(nod3!Eltseq(i)): i in RationalPoints(fl,QuadraticField(-3))] where nod3 is ChangeRing(nod,Qrm3);
// tangent at the node
tnode := Curve(P2,a-b);
// line through the flexes
l3 := Curve(P2,rep{x: x in DefiningPolynomials(fl)|Degree(x) eq 1});
pts := RationalPoints(fl,Qrm3);
pl1 := rep{x: x in Places(nod![0,1,0])|Degree(x) eq 1};
dm := DivisorMap(1*pl1);
_,md := IsInvertible(dm);
pts1 := [dm(p): p in pts];
P1 := Image(dm);
P1rm3 := ChangeRing(P1,Qrm3);
trs := [TranslationOfSimplex(P1rm3,[P1rm3(Qrm3)|Eltseq(pts1[i]): i in j]): j in [[1,2,3],[2,3,1]]]; 
tr3 := trs[1]^-1*trs[2];
nod3 := ChangeRing(nod,Qrm3);
dm3 := map<nod3->P1rm3|DefiningPolynomials(dm),DefiningPolynomials(md)>;
autnod3 := Extend(Normalization(Expand(dm3*tr3*dm3^-1)));
autnod3 := map<nod3->nod3|[nod3.1-rm3*nod3.2,-rm3*nod3.1+nod3.2,-9/4*nod3.1+3/4*rm3*nod3.2-2*nod3.3]>;

ex331 := Surface(P3w,w^2-(&*[CoordinateRing(P2)|DefiningPolynomial(i): i in [nod,l3,tnode,&join tans3]])@@to2);
// The surface is defined over Q, but the extra automorphism is only defined
// over Q(sqrt -3).  As far as I know the automorphism group of this surface
// over Q is generated by w -> -w.
ex331m3 := ChangeRing(ex331,Qrm3);
aut331m3 := map<e->e|Append([Evaluate(i/(1-rm3),[e.j: j in [1..3]]): i in DefiningPolynomials(autnod3)],e.4)> where e is ex331m3;
// an automorphism of order 3 ....
P3wm3 := AmbientSpace(ex331m3);
ls4 := LinearSystem(P3wm3,4);
act4 := [CoefficientMap(ls4)(i@@aut331m3): i in Sections(ls4)];
fix4 := [PolynomialMap(ls4)(i): i in Basis(Kernel(Matrix(act4)-1))];
// Since the automorphism and its conjugate are together defined over Q,
// the quotient is defined over Q, and it is a K3 surface.
P5 := ProjectiveSpace(Q,5);
to5 := map<ex331->P5|fix4>;
// in5 := Image(to5);
lss := [LinearSystem(P5,i): i in [2,3]];
repeat
  imi := to5(Scheme(ex331,[&+[Random(-3,3)*P3w.i: i in [1..3]]: j in [1,2]]));
  lss := [LinearSystem(l,imi): l in lss];
until [#Sections(i): i in lss] eq [2,15];
bs := Surface(P5,Ideal(&meet [BaseScheme(i): i in lss]));
to5 := map<ex331->bs|fix4>;
tok3 := MinimalModelKodairaDimensionZero(bs);
k3 := Codomain(tok3);

line := rep{x: x in PrimeComponents(Scheme(k3,k3.7))|Degree(x) eq 1};
P4 := ProjectiveSpace(Q,4);
to4 := map<k3->P4|DefiningPolynomials(line)>;
in4 := Image(to4);

// that was just to make point counting a bit faster
cts4 := [<p,p le 3 select 1 else #RationalPoints(ChangeRing(in4,GF(p)))-p^2-2*p-1>: p in PrimesUpTo(100)];
// this appears to be connected to 576.3.g.b, or certainly some twist of the
// form in level 16.
br331 := DefiningPolynomial(to2(Scheme(ex331,w)));
cts_big := [<p,p le 3 select 0 else &+[KroneckerSymbol(Evaluate(br331,Eltseq(i))): i in RationalPoints(ProjectiveSpace(GF(p),2))]>: p in PrimesUpTo(100)];
// the differences are 0 except for p congruent to 1 mod 12.
// There is some similarity to the coefficients of 576.3.g.f, but they are
// not the same.  The differences are
// [ <5, 0>, <13, 4>, <17, 0>, <29, 0>, <37, -52>, <41, 0>, <53, 0>, <61, -52>, <73, 92>, <89, 0>, <97, -4> ]
// and for p = 1 mod 12 we should have <13, -44>, <37, -52>, <61, -148>, <73, 92>, <97, -4>.
// however, no similarity past 100, so maybe just coincidence.

// example 3.3.2
cub := Curve(P2,b^2*c-a^3-c^3);  
dcub := Image(map<cub->P2|[Derivative(DefiningPolynomial(cub),i): i in [1..3]]>);
// and the reducible conic b^2 + 1/3 c^2 passes through 2 sets of 3 collinear
// cusps.
br332 := DefiningPolynomial(dcub)*(b^2+1/3*c^2);
ex332 := Surface(P3w,w^2-br332@@to2);

Pwnew := ProjectiveSpace(Q,[3,1,1,4]);
// the symmetry is obvious.
townew := map<ex332->Pwnew|[P3w.1^3,P3w.2,P3w.3,P3w.4]>;
// but in fact the quotient is rational.
cts332 := [<p,p le 3 select 1 else &+[KroneckerSymbol(Evaluate(br332,Eltseq(i))): i in RationalPoints(ProjectiveSpace(GF(p),2))]>: p in PrimesUpTo(100)];
// I don't understand these numbers.
