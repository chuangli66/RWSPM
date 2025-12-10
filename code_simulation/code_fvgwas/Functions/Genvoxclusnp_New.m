%%%%%%%%%%use top 200 SNP to appropriate top N0 SNP

function [pcorrect,ROI,rawpvalue]=Genvoxclusnp_New(N0,PX,X,YW,SNP,Wgmax,Ng,deta,index,sizeimg)

if nargin < 10
    Flag=0;
else
    Flag=1;
end

E1=eye(size(PX,1))-PX;
Q=SNP'*E1;%c*n
t1=sum(Q.^2,2);%c*1
t1=t1.^(-1);
t3=diag(t1);%c*c
t2=(Q*YW).^2;%c*V
deta=deta(:);
deta=deta';
t4=repmat(deta,[N0,1]);
W=t3*t2.*t4;%%%%%%C*V, real data's statistics

%%%%%%%%%%%%%%%association of voxels and SNPs
%%%%%%%%%%%real data statistic%%%%%
%%%%corrected pv
k1=mean(Wgmax);
k2=var(Wgmax);
k3=mean((Wgmax-k1).^3);
a=k3/(4*k2);
b=k1-2*k2^2/k3;
d=8*k2^3/k3^2;
pcorrect=1-cdf('chi2',(W-b)/a,d);%%%%%corrected pv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vv=size(PX,1)-size(X,2);
rawpvalue=1-cdf('f',W,1,vv);%%%%%raw pvalue


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Flag
    thred=0.005;
    alph=finv(1-thred,1,vv);
    rawpv=zeros(N0,length(index));
    rawpv(W(:)>alph)=1;
    ROI=cell(N0,1);
    for i=1:N0
        bwmask=zeros(sizeimg);
        bwmask(index)=rawpv(i,:);
        bwmask=bwlabeln(bwmask);
        STATS=regionprops(bwmask,'Area','Centroid');
        area=[STATS.Area];
        if ~isempty(area)
            for j=1:length(area)
                ROI{i}.pv(j)=length(find(Ng(:)>area(j)))/length(Ng);
            end
            ROI{i}.area=area';
            ROI{i}.bw=bwmask;
        end
    end
else 
    ROI=0;
end