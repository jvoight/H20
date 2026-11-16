// Count ALL lifts over congruence projective groups, including noncongruence lifts.
// Counts below are ACTUAL subgroups, not conjugacy classes.
SetColumns(0);
ChangeDirectory("sources");
load "csg1-lev48.dat";
ChangeDirectory("..");
PrintFile("all-lifts-over-congruence.txt","# CP g mu cusps c_reg all_pg2_lifts_per_fixed_projective_group projective_conjugates":Overwrite:=true);
audit:=AssociativeArray();
for line in Split(Read("audit.txt"),"\n") do
    fields:=Split(line," ");
    if #fields ne 4 or fields[1] eq "#" then continue; end if;
    audit[fields[1]]:=StringToInteger(fields[3]);
end for;
total:=0; cong_total:=0;
for r in L do
    if r`index gt 36 or 1 in r`c2 then continue; end if;
    e3:=#[x:x in r`c3|x eq 1];
    cr:=(r`index+2*e3-12)/3; c:=#r`cusps;
    if cr lt 0 or cr gt c or Denominator(cr) ne 1 then continue; end if;
    cr:=Integers()!cr;
    assert IsEven(cr);
    // A projective orbifold of signature (g;3^e3;c) has
    // 2^(2g) lifts for each admissible cusp-sign pattern;
    // precisely the patterns with an even number of regular cusps occur.
    n:=2^(2*r`genus)*Binomial(c,cr);
    copies:=r`GLZConj*r`length;
    total+:=n*copies;
    assert IsDefined(audit,r`name) and audit[r`name] le n;
    cong_total+:=audit[r`name]*copies;
    PrintFile("all-lifts-over-congruence.txt",Sprintf("%o %o %o %o %o %o %o",r`name,r`genus,r`index,c,cr,n,copies));
    if r`name eq "4G0" then
        assert n eq 15 and copies eq 1 and audit[r`name] eq 3;
        print "Above projective Gamma(4):",n,"lifts with pg=2";
    end if;
end for;
assert cong_total eq 558 and total eq 1745;
print "All actual pg=2 subgroups with congruence projective image:",total;
print "Congruence actual subgroups (separate complete enumeration):",558;
print "Noncongruence actual subgroups with congruence projective image:",total-558;
print "PASS: all-lift count from the orbifold presentation";
quit;
