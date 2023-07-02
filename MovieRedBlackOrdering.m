clc;
clear all;
close all;

%Movie for red-black ordering algorithm
N=3;
A=sparse(toeplitz([2,-1,zeros(1,N-2)]));
I=speye(N,N);
B=kron(A,I)+kron(I,A);
[x,y]=ndgrid(1:N,1:N);
p=redblack(B);
B=B(p,p);
R=chol(B(p,p));
xy=[x(p);y(p)]';
n=size(B,1);
L=zeros(n);
lpars={'marker','.','linestyle','none','markersize',64};
tpars={'fontname','helvetica','fontsize',16,'horiz','center','col','g'};
for i=1:n
 clf
 gplot(B,xy);
 line(xy(i:N^2,1),xy(i:N^2,2),'color','k',lpars{:});
 for j=i:n
 degree=length(find(B(:,j)))-1;
 text(xy(j,1),xy(j,2),int2str(degree),tpars{:});
 end
 axis equal,axis off,axis([1,N,1,N])
 pause(0.3);
 line(xy(i,1),xy(i,2),'color','r',lpars{:});
 pause(0.3);
 L(i,i)=sqrt(B(i,i));
 L(i+1:n,i)=B(i+1:n,i)/L(i,i);
 B(i+1:n,i+1:n)=B(i+1:n,i+1:n)-L(i+1:n,i)*L(i+1:n,i)';
 B(i, i:n)=0;
 B(i:n, i)=0;
end
spy(L,32) 

function p=redblack(A)
%REDBLACK: Computes the permutation vector implementing the red-black
%ordering algorithm 
n=size(A,1);
temp = 1:n;
count = 0;
odd = 2;
flag =1;
for m = 1:n
 if flag
 p(m)=temp(m+count);
 count = count+1;
 if p(m)==n|p(m)+2>n
 flag = 0;
 end
 else
 p(m)=temp(odd);
 odd = odd+2;
 end 
end
end