%% function FVGWAS_GSIS_server
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Matrix prepraring and Data preprocessing



% s is  the number of Setting 1-7
for s = 1:7
    % File loading path
    root = "/.../Setting" + num2str(s);
    % file storage path
    save_root = "/.../Setting"+ num2str(s) + "/";
    n=200;
    for l = 1:100
        snp_path =  root + '/sim_' +  num2str(l) +'/x.csv';
    
        SNP1 = readtable(snp_path);
        SNP1 = SNP1(2:end,2:end);
        SNP1 = table2array(SNP1);
        SNP1 = [SNP1,randi([0, 2], n, 100)];
        X = randn(n,5);
        p_vec = [];
        img_path = root +  '/sim_'+ num2str(l) + '/Y_slide'+ num2str(l) + '.mat';
        for region_num = 1:6534
            load(img_path);
            YW = reshape( Y(region_num,:),200,9);
    
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
            p_vec(region_num) = pv(1);
        end
        l
        ind_name = strcat(save_root,'/test_fvgwas',num2str(l),'.csv');
        writematrix([(1:6534)',p_vec'],ind_name);
    end
end