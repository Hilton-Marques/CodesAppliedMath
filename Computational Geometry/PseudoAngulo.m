
v = [4,-2];
out = getPseudoAngle(v)

function out = getPseudoAngle(v)
if v(1) > 0
    if v(2) > 0
        if v(2) > v(1)
            out = (2 - v(1)/v(2));
            return
        end
        out = v(2)/v(1);
    else
        if -v(2) > v(1)
            out = 6  + v(1)/-v(2);
            return
        end
        out = 8 - (-v(2)/v(1));
    end
else
    if v(2) > 0
        if v(2) > -v(1)
            out = 2 + -v(1)/v(2);
            return
        end
        out = 4 - (v(2)/-v(1));
    else
        if -v(2) > -v(1)
            out = 6 - (-v(1)/-v(2));
            return
        end
        out = 4 + (-v(2)/-v(1));
    end
end
end

