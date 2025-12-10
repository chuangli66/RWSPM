function [pp,pv, W]=GenBraAsso_New(SNP,YW,PX,deta)


E3=eye(size(PX,1))-PX;
Q=SNP'*E3;
t1=sum(Q.^2,2);
t1=t1.^(-1);
E4=YW*(repmat(deta,[1,size(YW,1)]).*YW');
[U,D]=eig(E4);
D = real(D);
D(D(:)<0)=0;
t4=U*(D^(1/2))*U';
t4=real(t4);
t2=sum((Q*t4).^2,2);%%%%%%C*1
W=t1.*t2;%%%%%%%
W=W./length(deta);

k1=mean(W);
k2=var(W);
k3=mean((W-k1).^3);
a=k3/(4*k2);
b=k1-2*k2^2/k3;
d=8*k2^3/k3^2;
pv=1-cdf('chi2',(W-b)/a,d);
pp=-log10(pv);