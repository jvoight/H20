// Independent standard-family checks against the saved representatives.
load "representatives.m";
function MP(m,a)
 return Sym(m^2)![1+(((i mod m)*a[1]+(i div m)*a[3]) mod m)+m*(((i mod m)*a[2]+(i div m)*a[4]) mod m):i in [0..m^2-1]];
end function;
function Ent(m,p)
 u:=2^p-1;v:=(m+1)^p-1;return [u mod m,u div m,v mod m,v div m];
end function;
for ex in [<12,-4>,<14,-7>,<15,-3>,<15,-15>,<16,-8>,<18,-3>,<9,0>,<10,0>] do
 N,D:=Explode(ex);
 G:=sub<Sym(N^2)|MP(N,[0,1,-1,0]),MP(N,[1,1,0,1])>;
 els:=[];
 for x in G do
  a:=Ent(N,x);
  if a[3] eq 0 and ((D eq 0 and a[1] eq 1 and a[4] eq 1) or (D ne 0 and KroneckerSymbol(D,a[4]) eq 1)) then Append(~els,x); end if;
 end for;
 H:=sub<G|els>; hits:=[];
 for i in [1..#CongruencePG2] do
  r:=CongruencePG2[i];
  if r[4] ne N then continue; end if;
  K:=sub<G|[MP(N,a):a in r[11]]>;
  if IsConjugate(G,H,K) then Append(~hits,<i,r[1]>); end if;
 end for;
 assert #hits eq 1;
 print ex,hits;
end for;
print "PASS: all six quadratic refinements and Gamma_1(9), Gamma_1(10)";
quit;
