classdef SphericalHarmonics < handle
    properties
        n
        m
        MValuesPos
        values
    end
    methods
        function this = SphericalHarmonics(n)
            if nargin > 0
                this.n = n;
                this.m = 2*(n-1)+1;
                mPos = (this.m - 1)/2 + 1;
                this.MValuesPos = zeros(1,mPos);
            end
        end
        function addNNValue(this,value)
            this.MValuesPos(end) = value;
        end
        function addMValues(this,value,j)
            this.MValuesPos(j) = value;
        end
        function merge(this)
            mPos = (this.m - 1)/2 + 1;
            this.values = zeros(1,mPos-1);
            for i = mPos:-1:2
                value = this.MValuesPos(i);
                this.values(mPos - i + 1) = ((-1)^(i-1))*conj(value);
            end
            this.values = [this.values,this.MValuesPos];
        end
    end
end