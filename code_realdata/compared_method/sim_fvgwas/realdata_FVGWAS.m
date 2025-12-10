%% function FVGWAS_GSIS_server
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = readtable('E:/Region_based_function/zhu_code/P1CodeDocs/P1CodeDocs/realdata/biom12748-sup-0002-suppinfo-s2/hippocampussurface.dat');

%% load the snp data SNP
snp_path = '/.../ADNI_1_3.csv';
SNP = readtable(snp_path);
SNP1 = table2array(SNP(:,2:end));
%% load the covariate data X
cov_inf_path ='/.../covariates_inf.csv';
X = readtable(cov_inf_path);
X = table2array(X(:,2:end));
%% load the image data YW
% Y1= load('/data/left_res_ori.mat');
Y1= load('/data/right_res_ori.mat');
YW = Y1.Y;
%%%%%%%%%%%%%%%%
lamda=0.001;%%%%% parameters to avoid matrix sigularity
PX=X/(X'*X+lamda*eye(size(X,2)))*X';
n=size(X,1);%%%%%number of individual
p=size(X,2);%%%%feature dimention of X

%%%%%%%%% Preprocess SNP
num=zeros(3,size(SNP1,2));
for i=1:3
    bw=zeros(size(SNP1));
    bw(SNP1(:)==i-1)=1;
    num(i,:)=sum(bw);
end
[~,maxnum]=max(num);
maxnum=maxnum-1;
[r,c]=find(SNP1<0);
for i=1:length(r)
    SNP1(r(i),c(i))=maxnum(c(i));
end

%%%%%%%%%%  Calculate Sigma_hat across all loci  %%%%%%%%%
BatchSize=5000;
deta=zeros(size(YW,2),1);
if size(YW,2)>BatchSize
    for kk=0:fix(size(YW,2)/BatchSize)
        if kk<fix(size(YW,2)/BatchSize)
            temp=YW(:,BatchSize*kk+1:BatchSize*(kk+1))'*(eye(size(PX,1))-PX)*YW(:,BatchSize*kk+1:BatchSize*(kk+1))/(n-p);
            deta(BatchSize*kk+1:BatchSize*(kk+1))=diag(temp);
        else
            if BatchSize*kk+1>size(YW,2)
                break;
            end
            temp=YW(:,BatchSize*kk+1:end)'*(eye(size(PX,1))-PX)*YW(:,BatchSize*kk+1:end)/(n-p);
            deta(BatchSize*kk+1:end)=diag(temp);
        end
    end
else
    temp=YW'*(eye(size(PX,1))-PX)*YW/(n-p);
    deta=diag(temp);
end
deta(deta(:)~=0)=deta(deta(:)~=0).^(-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pp,pv]=GenBraAsso_New(SNP1,YW,PX,deta);
pv = [right_top_index,pv];
% ind_name = strcat('/.../fvgwas_global_left.txt');
ind_name = strcat('/.../fvgwas_global_right.txt');
save(ind_name,'pv','-ascii');