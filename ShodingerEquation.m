harmonic(10,10,8)
function harmonic(L,n,k)
h = 1/n; N =2*n*L+1;
K = toeplitz([2 -1 zeros(1,N-2)]);
H = K/(h^2) + diag((-L:h:L).^2);
[V,F] = eig(H);
E = diag(F);
E = E(1:k);
j = 1:k; plot(j,E)
end
