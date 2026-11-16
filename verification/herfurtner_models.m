// Reconstructed from rows 7, 13, 15, 25, 28 of legacy-surfaces/maximal.m.
// Return the untwisted equation y^2 = x^3 - g2*x + 2*g3 and its twist.
// The twist is the product of the monic finite bad-place polynomials.
// Its constant squareclass is part of the arithmetic model: do not rescale it.
function HerfurtnerInput(id)
    Q<t> := FunctionField(Rationals());
    case id:
    when 7:
        g2 := 3*(8*t-1)*(8*t^3+87*t^2+96*t-64);
        g3 := (8*t-1)*(64*t^5+1040*t^4+3925*t^3+100*t^2+3200*t-512);
        d := t^3+79/8*t^2-5/4*t;
    when 13:
        g2 := 3*t^2*(9*t^2-8);
        g3 := t^2*(27*t^4-36*t^2+8);
        d := t^3-t;
    when 15:
        g2 := 3*(t-1)^2*(9*t^2+14*t+9);
        g3 := (t-1)^2*(27*t^4+36*t^3+2*t^2+36*t+27);
        d := t^3-t;
    when 25:
        g2 := 3*t*(t^2+6*t-3);
        g3 := 6*t^2*(t^2+3);
        d := t^3-6*t^2-3*t;
    when 28:
        g2 := 3*(t-1)^2*(t+1)^2;
        g3 := (t-1)^2*(t+1)^2*(t^2+1);
        d := t^3-t;
    else
        error "Supported Herfurtner indices are 7, 13, 15, 25, 28";
    end case;
    return EllipticCurve([-g2,2*g3]), d;
end function;
