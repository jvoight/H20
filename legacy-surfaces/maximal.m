// we can get elliptic surfaces of Picard number 30 from quadratic twists
// of the surfaces in Herfurtner that have no stars.

// The D8D8D6D6 is of course the quadratic twist of the Legendre family.
// The L-function factors as p_2^2 for p = 1 mod 4 and is irreducible for
// p = 3 mod 4.  Similarly there is a D12D6D5D5 which is 2-isogenous to that.

// The D12E6D5D5 has L-function factoring as p_2 p_1 p_1' for p = 5, 7 mod 12
// and is irreducible for p = 1, 11 mod 12.

Qs<s> := FunctionField(Rationals());
gis := [[12*(s^4-s^2+1),4*(2*s^6-3*s^4-3*s^2+2)], // Legendre
        [3*s*(9*s^3-8),27*s^6-36*s^3+8], // I9I1I1I1
        [12*(s^4-4*s^3+2*s+1),4*(2*s^6-12*s^5+12*s^4+14*s^3+3*s^2+6*s+2)], // I6I3I2I1
        [3*(s^4-12*s^3+14*s^2+12*s+1),s^6-18*s^5+75*s^4+75*s^2+18*s+1], // I5I5I1I1]
        // I leave out the isogenous ones, namely I3I3I3I3 and I8I2I1I1
        [12*s*(s^3-6*s^2+15*s-12),4*s*(2*s^5-18*s^4+72*s^3-144*s^2+135*s-27)], // II I8I1I1
        [12*s*(9*s^3+36*s^2+42*s+14),12*s*(18*s^5+108*s^4+234*s^3+222*s^2+87*s+8)], // II I7I2I1
        [3*(8*s-1)*(8*s^3+87*s^2+96*s-64),(8*s-1)*(64*s^5+16*5*13*s^4+25*157*s^3+100*s^2+2^7*5^2*s-2^9)], // II I5I4I1
        [3*(s-3)*(81*s^3-9*s^2-53*s-27),(s-3)*(3^6*s^5-3^5*5*s^4-2*3^3*5^2*s^3-350*s^2-3^3*5^2*s-243)], // II I5I3I2
        [12*s*(s^3+4*s^2+10*s+6),4*s^2*(2*s^4+12*s^3+42*s^2+70*s+63)], // III I7 I1I1
        [12*s*(s^3-6*s^2+9*s-3),4*s^2*(2*s^4-18*s^3+54*s^2-63*s+27)], // III I6I2I1
        [75*(5*s-1)*(5*s^3+45*s^2+39*s-25),25*(5*s-1)^2*(25*s^4+340*s^3+2*3*181*s^2+100*s+5^4)], // III I5I3I1: typo in Herfurtner
        [3*(s-1)*(16*s^3-3*s-1),(s-1)^2*(64*s^4+32*s^3+6*s^2+5*s+1)], // III I4I3I2
        [3*s^2*(9*s^2-8),s^2*(27*s^4-36*s^2+8)], // IV I6I1I1
        [12*s^2*(s^2+8*s+10),4*s^2*(2*s^4+24*s^3+78*s^2+56*s+27)], // IV I5I2I1
        [3*(s-1)^2*(9*s^2+14*s+9),(s-1)^2*(27*s^4+36*s^3+2*s^2+36*s+27)], // IV I3I3I2
        [3*(s^2-13*s+49)*(s^2-5*s+1),(s^2-13*s+49)*(s^4-14*s^3+63*s^2-70*s-7)], // II II I7I1
        [3*(s-1)*(s+1)*(9*s^2-1),(s-1)*(s+1)*(27*s^4-18*s^2-1)], // II II I6I2
        [12*s*(2*s-1)*(s+4),2*(2*s-1)*(s+4)*(s^4+4*s^3+8*s-4)], // II II I4I4
        // [2*s*(s-1)*(6*s^2+6*z*s-(3+z)),2*s^2*(s-1)*(4*s^3-2*(1-3*s)*s^2-4*(2+z)*s+(5-z))], // III II I6I1 where z = sqrt(-3)
        [3*s*(2*s-27)*(2*s^2-35*s+140),s*(2*s-27)^2*(2*s^3-39*s^2+222*s-250)], // III II I5I2
        [3*(s-1)*(s+27)*(16*s^2+80*s-243),(s-1)^2*(s+27)*(64*s^3+32*43*s^2+2*3^5*s+3^9)], // III II I4I3
        [3*s^2*(s-3)*(s+5),s^2*(s-3)*(s^3+6*s^2-3*s-32)], // IV II I5I1
        [36*s^2*(s-1)*(3*s-1),4*s^2*(s-1)*(54*s^3-54*s^2+9*s-1)], // IV II I4I2
        [3*(s^2+1)*(s^2+6*s+4),(s^2+1)^2*(s^2+9*s+19)], // III III I5I1
        [3*(s-1)*(s+1)*(4*s^2-1),(s-1)^2*(s+1)^2*(8*s^2+1)], // III III I4I2
        [3*s*(s^2+6*s-3),6*s^2*(s^2+3)], // III III I3I3
        [12*s^2*(s-1)*(s+5),4*s^2*(s-1)^2*(2*s^2+16*s+27)], // IV III I4I1
        [3*s^2*(s-1)*(9*s-5),s^2*(s-1)^2*(27*s^2-9*s+2)], // IV III I3I2
        [3*(s-1)^2*(s+1)^2,(s-1)^2*(s+1)^2*(s^2+1)], // IV IV I2I2
        [3*s^2*(s-1)*(s+1),s^2*(s-1)^2*(s+1)^2], // IV III III I2
        [3*s*(s^3-1),(s^3-1)^2], // III III III I3
        [3*s^2*(s-1)*(s+3),s^2*(s-1)^2*(s+3)*(s+2)], // IV III II I3
        [12*s^2*(s-1)*(s+1),4*s^2*(s-1)*(s+1)*(2*s^2-1)], // IV II II I4
        [3*(s^2+2)*(s+5)*(s+1),(s^2+2)^2*(s+5)*(s+4)], // III III II I4
        [3*s*(s^2+11*s+64)*(s+3),s^2*(s^2+11*s+64)*(s^2+10*s+45)], // III II II I5
        [12*s*(s^3-1),4*(s^3-1)*(2*s^3-1)]]; // II II II I6
ecs := [EllipticCurve([-x[1],2*x[2]]): x in gis];
tws := [QuadraticTwist(e,&*[Numerator(Qs!Minimum(i[1])): i in LocalInformation(e)]): e in ecs];
lfss := [];
for e in tws do
  lfs := [];
  for p in PrimesUpTo(200) do try
      ep := ChangeRing(e,FunctionField(GF(p)));
      lf := LFunction(ep);
      Append(~lfs,<p,lf>);
      printf "."; 
    catch err; end try;
  end for;
  " done curve with Kodaira symbols", KodairaSymbols(e);
  Append(~lfss,lfs);
end for;

// Q(i): 1, 23, 28, 32
// 1 and 23: https://www.lmfdb.org/L/4/2e10/1.1/c2e2/0/0
// 28: https://www.lmfdb.org/L/4/288e2/1.1/c2e2/0/1
// 32: ? (can't find anything suitable with a37=76)
// anyway, something weird is going on--one kind of L-factor for
// 5, 17, 29, 37, 41 and another for the other primes?

// Q(sqrt -3): 2, 3
// 2: ? though Euler factors look right
// 3: http://www.lmfdb.org/L/4/24e2/1.1/c2e2/0/0

// 4: ?
// 5: ?
// 6: ? (5 and 6 both have factors of p^2T^2-1 in the L-function for
// many but not all primes)
// 13: again, one kind of L-factor at 5, 17, 29, 37, 41 (and 53 61 89;
// what's wrong with 13 73 97?  They are 1 mod 12, but so is 37.)
// 14: I have no idea what's going on, but it looks like the L-function is
// a function of T^2 if and only if p is not 3, 7, 9 mod 20
// 17: http://www.lmfdb.org/L/4/640e2/1.1/c2e2/0/11 (quite different-looking
// L-series)


// so I have explained 5 of them, namely:
// d6d6d8d8 (1), d9e7e7d5 (23)
// e8e8d6d6 (28)
// d10d7d6d5 (3)
// e6e6d10d6 (17)

// Looks like 10 and 12, 13 and 15 are quadratic twists of each other
// where 10, 12 is explained by 2-isogeny and 13, 15 by 3-isogeny.

// note also: 11, 26, 27 are bad at 5 and 19, 20 are bad at 7
// everything is bad at 2, 3, and that's it.
