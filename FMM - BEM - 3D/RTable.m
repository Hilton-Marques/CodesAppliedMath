classdef RTable < handle
    properties
        % A table with the R's functions evaluated in a y vector
        n
        Ry
    end
    methods
        function this = RTable(n,y)
            this.n = n;
            yb = y(1) + 1i*y(2);
            conjY = conj(yb);
            AbsSquaredY = yb*conjY;
            y3 = y(3);
            if n < 7
                Ry = [1,(1/2)*conjY,y3,(-1/2)*yb,(1/8)*conjY^2,(1/2)*conjY*y3,( ...
                    1/4)*((-1)*AbsSquaredY+2*y3^2),(-1/2)*y3*yb,(1/8)*yb^2,( ...
                    1/48)*conjY^3,(1/8)*conjY^2*y3,(-1/16)*conjY*(AbsSquaredY+( ...
                    -4)*y3^2),(-1/4)*AbsSquaredY*y3+(1/6)*y3^3,(1/16)*( ...
                    AbsSquaredY+(-4)*y3^2)*yb,(1/8)*y3*yb^2,(-1/48)*yb^3,( ...
                    1/384)*conjY^4,(1/48)*conjY^3*y3,(-1/96)*conjY^2*( ...
                    AbsSquaredY+(-6)*y3^2),(1/48)*conjY*y3*((-3)*AbsSquaredY+4* ...
                    y3^2),(1/192)*(3*AbsSquaredY^2+8*(y3^4+(-3)*y3^2*yb^2)), ...
                    (-1/48)*y3*((-3)*AbsSquaredY+4*y3^2)*yb,(-1/96)*( ...
                    AbsSquaredY+(-6)*y3^2)*yb^2,(-1/48)*y3*yb^3,(1/384)*yb^4, ...
                    (1/3840)*conjY^5,(1/384)*conjY^4*y3,(-1/768)*conjY^3*( ...
                    AbsSquaredY+(-8)*y3^2),(1/96)*conjY^2*y3*((-1)*AbsSquaredY+ ...
                    2*y3^2),(1/384)*(AbsSquaredY^2+(-12)*AbsSquaredY*y3^2+8* ...
                    y3^4)*yb,(1/960)*(15*AbsSquaredY^2*y3+8*(y3^5+(-5)* ...
                    y3^3*yb^2)),(-1/384)*(AbsSquaredY^2+(-12)*AbsSquaredY* ...
                    y3^2+8*y3^4)*yb,(1/96)*y3*((-1)*AbsSquaredY+2*y3^2)* ...
                    yb^2,(1/768)*(AbsSquaredY+(-8)*y3^2)*yb^3,(1/384)*y3* ...
                    yb^4,(-1/3840)*yb^5,(1/46080)*conjY^6,(1/3840)*conjY^5*y3, ...
                    (-1/7680)*conjY^4*(AbsSquaredY+(-10)*y3^2),(1/2304)* ...
                    conjY^3*y3*((-3)*AbsSquaredY+8*y3^2),(1/3072)*AbsSquaredY* ...
                    (AbsSquaredY^2+(-16)*AbsSquaredY*y3^2+16*y3^4),(1/1920)* ...
                    y3*(5*AbsSquaredY^2+(-20)*AbsSquaredY*y3^2+8*y3^4)*yb,( ...
                    1/11520)*((-5)*AbsSquaredY^3+90*AbsSquaredY^2*y3^2+(-120)* ...
                    AbsSquaredY*y3^4+16*y3^6),(-1/1920)*y3*(5*AbsSquaredY^2+( ...
                    -20)*AbsSquaredY*y3^2+8*y3^4)*yb,(1/3072)*yb^2*( ...
                    AbsSquaredY^2+16*(y3^4+(-1)*y3^2*yb^2)),(-1/2304)*y3*(( ...
                    -3)*AbsSquaredY+8*y3^2)*yb^3,(-1/7680)*(AbsSquaredY+(-10)* ...
                    y3^2)*yb^4,(-1/3840)*y3*yb^5,(1/46080)*yb^6];
            end
            this.Ry = Ry;
        end
        function out = pos(this,position)
            if (position <= 0)
                out = 0;
                return;
            end
            out = this.Ry(position);
        end
        function out = getRow(this,n)
            posAfter = n^2;
            out = this.Ry(1:posAfter);
        end
    end
end