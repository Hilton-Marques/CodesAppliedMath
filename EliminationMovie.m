clc;
clear all;
close all;
N=3;
A=sparse(toeplitz([2,-1,zeros(1,N-2)]));
I=speye(N,N);
B=kron(A,I)+kron(I,A);
[x,y]=ndgrid(1:N,1:N);
B0=B;
L0=chol(B0)';

% p=symamd(B);
p=realmmd(B);
B=B(p,p);
R=chol(B);
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
  pause(.01);

  line(xy(i,1),xy(i,2),'color','r',lpars{:});
  pause(.01);

  L(i,i)=sqrt(B(i,i));
  L(i+1:n,i)=B(i+1:n,i)/L(i,i);
  B(i+1:n,i+1:n)=B(i+1:n,i+1:n)-L(i+1:n,i)*L(i+1:n,i)';
  B(i, i:n)=0;
  B(i:n, i)=0;
end
figure,spy(B0)
figure,spy(B)
figure,spy(L0,32)
figure,spy(L,32)