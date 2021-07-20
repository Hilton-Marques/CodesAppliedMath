key = '0b010001000111';

level = numel(key)/3;
keyXYZ = decode(key);
vicinityBase = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
    (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1)-1),(keyXYZ(1));...
    (keyXYZ(2)-1),(keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),...
    (keyXYZ(2)+1),(keyXYZ(2)+1),(keyXYZ(2)+1),keyXYZ(2),keyXYZ(2);
    (keyXYZ(3)-1), (keyXYZ(3)-1), (keyXYZ(3)-1), (keyXYZ(3)-1),(keyXYZ(3)-1), ...
    (keyXYZ(3)-1),(keyXYZ(3)-1),(keyXYZ(3)-1),(keyXYZ(3)-1)];
vicinityMiddle = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
    (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1)-1);...
    (keyXYZ(2)-1),(keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),...
    (keyXYZ(2)+1),(keyXYZ(2)+1),(keyXYZ(2)+1),keyXYZ(2);
    (keyXYZ(3)), (keyXYZ(3)), (keyXYZ(3)), (keyXYZ(3)),(keyXYZ(3)), ...
    (keyXYZ(3)),(keyXYZ(3)),(keyXYZ(3))];
vicinitySuper = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
    (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1)-1),(keyXYZ(1));...
    (keyXYZ(2)-1),(keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),...
    (keyXYZ(2)+1),(keyXYZ(2)+1),(keyXYZ(2)+1),keyXYZ(2),keyXYZ(2);
    (keyXYZ(3)+1), (keyXYZ(3)+1), (keyXYZ(3)+1), (keyXYZ(3)+1),(keyXYZ(3)+1), ...
    (keyXYZ(3)+1),(keyXYZ(3)+1),(keyXYZ(3)+1),(keyXYZ(3)+1)];
bin = encode([vicinityBase,vicinityMiddle,vicinitySuper]);
check = decode(bin);

function out = decode(key)
binx = strcat('0b',key(:,5:3:end));
biny = strcat('0b',key(:,4:3:end));
binz = strcat('0b',key(:,3:3:end));
out = [int16(str2num(binx)),int16(str2num(biny)),int16(str2num(binz))];
end
function out = vicinityIndex(key)
keyDec = base2dec(key,4);
keyXYZ = this.decode(keyDec);
vicinity = [(keyXYZ(1)-1),keyXYZ(1),(keyXYZ(1)+1),(keyXYZ(1)+1),...
    (keyXYZ(1)+1),keyXYZ(1),(keyXYZ(1)-1),(keyXYZ(1) -1);(keyXYZ(2)-1),...
    (keyXYZ(2)-1),(keyXYZ(2)-1),keyXYZ(2),(keyXYZ(2)+1),(keyXYZ(2)+1),...
    (keyXYZ(2)+1),keyXYZ(2)];
out = this.encode(vicinity');
end
function out = encode(bin)
x = bin(1,:);
y = bin(2,:);
z = bin(3,:);
xBin = dec2bin(x);
yBin = dec2bin(y);
zBin = dec2bin(z);
nEl = max([numel(xBin(1,:)),numel(yBin(1,:)),numel(zBin(1,:))]);
xBin = dec2bin(x,nEl);
yBin = dec2bin(y,nEl);
zBin = dec2bin(z,nEl);
out = char(zeros(26,3*nEl));
out(:,3:3:end) = xBin;
out(:,2:3:end) = yBin;
out(:,1:3:end) = zBin;
out = strcat('0b',out);
end