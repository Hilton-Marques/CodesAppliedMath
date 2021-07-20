% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
points = [1,0,0;...
    1,1,0;...
    0,1,0;...
    0,0,0;...
    1,0,1;...
    1,1,1;...
    0,1,1;...
    0,0,1];
Faces(6) = Face();
Faces(1) = Face(points([1,2,6,5],:),1,'0');
Faces(2) = Face(points([2,3,7,6],:),1,'0');
Faces(3) = Face(points([3,4,8,7],:),1,'0');
Faces(4) = Face(points([1,5,8,4],:),1,'0');
Faces(5) = Face(points([1,4,3,2],:),1,'0');
Faces(6) = Face(points([5,6,7,8],:),1,'0');
lmax = 3;
Campo = Field(lmax,Faces);
elementos = Campo.cells;
hash = Campo.levels;
level1 = elementos(logical(hash(:,1)));
level2 = elementos(logical(hash(:,2)));
level3 = elementos(logical(hash(:,lmax)));

hold on
view(45,45);
for i = 1:3
    cell = level1{i};
    cell.show(1);
    for i = 1:4
        child = cell.child{i};
        show(child);
    end
end
keyboard
function show(cell)
cell.show(cell.level*0.5);
if ~cell.isDivided
    return
end
for i = 1:4
    child = cell.child{i};
    show(child);
end
end



