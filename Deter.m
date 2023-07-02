x1 =[1 2];
x2 = [3 4];
x3 = [3 5];

A = [x1 1;x2 1;x3 1]
A1 = [x1;x2];
A2 = [x1;x3];
A3 = [x3;x1];
det(A)
det(A1)
det(A2)
det(A3)
B = [(x2-x1)' (x3-x1)']
det(B)
