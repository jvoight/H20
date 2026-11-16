h02 := function(d,n);
  series := 1/((1-z)*(1+z*y))*&*[((1+z*y)^d[i] - (1-z)^d[i])/((1+z*y)^d[i]+y*(\
1-z)^d[i]) : i in [1..#d]];
  chi0 := Coefficient(series,n+#d);
  return (-1)^n*Coefficient(chi0,0) + 1;
end function;
 