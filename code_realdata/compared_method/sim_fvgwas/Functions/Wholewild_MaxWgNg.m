
function [Wgmax,Ngmax]=Wholewild_MaxWgNg(N0,X,PX,YW,deta,ZZ,G,index,sizeimg)

if nargin < 9
    Flag=0;
else
    Flag=1;
end

E1=eye(size(PX,1))-PX;
QQ=ZZ'*E1;
tt1=sum(QQ.^2,2);
tt1=tt1.^(-1);

ebeta=(X'*X)\X'*YW;
eerorr=YW-X*ebeta;
deta=deta(:);
E4=eerorr*(repmat(deta,[1,size(YW,1)]).*eerorr');
[U,D]=eig(E4);
D = real(D);
D(D(:)<0)=0;
tt4=U*(D^(1/2))*U';
tt4=real(tt4);

P=eerorr.^2;%%%%n*V

% Wg=zeros(size(tt1,1),G);
Wgmax=zeros(G,1);
roinum=1;
Ngmax=zeros(G*roinum,1);
thred=0.005;
vv=size(PX,1)-size(X,2);
alph=finv(1-thred,1,vv);
% step=1;
rand('seed',sum(100*clock));
% k=0;
for i=1:G
%     tic
    eta=(rand(size(eerorr,1),1)>0.5)*2-1;
% %     eta=normrnd(0,1,[size(eerorr,1),1]);
    eta1=diag(eta);
    t2=sum((QQ*eta1*tt4).^2,2);%%%%C*1
    Wg=tt1.*t2./length(deta);%%%%C*1
    [~,indx]=sort(Wg,'descend');
    SNP=ZZ(:,indx(1:N0));
       
    Q=SNP'*E1;
    temp=Q.^2;
    t1=temp*P;%%%%%%c*V
    t1=t1.^(-1);
    t2=Q*(eerorr.*repmat(eta,[1,size(eerorr,2)]));
    temp=t1.*(t2.^2);%%%%N0*V
    Wgmax(i)=max(temp(:));

    if Flag
        rawpv=zeros(N0,length(index));
        rawpv(temp(:)>alph)=1;
        
        clear t2 temp t1 Q SNP
        Nc=zeros(N0,1);
        for jj=1:N0
            bwmask=zeros(sizeimg);
            bwmask(index)=rawpv(jj,:);
            bwmask=bwlabeln(bwmask);
            STATS=regionprops(bwmask,'Area');
            area=[STATS.Area];
            if ~isempty(area)
                Nc(jj)=max(area);
            end
            clear bwmask area STATS 
        end
        
        Ngmax(i)=max(Nc);

    end
    i
%     toc
end
