function [fM,fZe,fphi,LM] = Sim_CGM(fM,max_iter,fpop,A,ftM,fZ)
    N = size(fM,1);
    alpha = 2*ones(1,1,1,size(fM,5));
    theta = ones(size(fM,2),size(fM,3),size(fM,4),size(fM,5));
    theta = normalize(theta,A);
    %Construct initial phi
    fphi = zeros(size(fM));
    fmf = ones(N,size(fM,2),size(fM,3),size(fM,4));
    fmb = fmf;
    for i=1:N
        M = squeeze(fM(i,:,:,:,:));
        Z = squeeze(fZ(i,:,:,:));
        mf = squeeze(fmf(i,:,:,:));
        mb = squeeze(fmb(i,:,:,:));
        fphi(i,:,:,:,:) = init_phi(M,A);
        fM(i,:,:,:,:) = update_m(M,Z,mf,mb,squeeze(fphi(i,:,:,:,:)),A);
    end
    conv = 0;
    %Loop till convergence 
    iter = 0;
    LM=zeros(N,max_iter,3);
    prob_diff=[];
    beta = 0.9;
    while conv == 0
        iter = iter+1;
        %For each dataset...
        for f=1:N
            mf = squeeze(fmf(f,:,:,:));
            mb = squeeze(fmb(f,:,:,:));
            M = squeeze(fM(f,:,:,:,:));
            Z = squeeze(fZ(f,:,:,:));
            pop = squeeze(fpop(f));
            tM = squeeze(ftM(f,:,:,:,:));
            
            C = const_mult(M,Z,A,pop);
            phi = theta.*C;
            
            %Calculate messages
            [mf,mb] = calc_messages(mf,mb,phi,A);
        
            %Update M
            M_n = update_m(M,Z,mf,mb,phi,A);
            M_n = beta*M + (1-beta)*M_n;


            Ze = mf.*mb;
            for t=1:size(M,3)
                Ze(:,:,t) = pop*Ze(:,:,t)/sum(sum(Ze(:,:,t)));
            end
 
            M = M_n;
            acv = Mdiff(tM,M);
            acc(f,:) = acv;
            cct = check_constraints(Ze,M,A);
            LM(f,iter,:) = [acv(1),pop,cct];
            
            fphi(f,:,:,:,:) = phi;
            fmf(f,:,:,:) = mf;
            fmb(f,:,:,:) = mb;
            fM(f,:,:,:,:) = M;
            fZe(f,:,:,:) = Ze;
        end
        to = theta;
        ao = alpha;

        %Adjust alpha and theta values
        theta = alpha+squeeze(sum(fM,1))+.01;
        theta = exp(psi(theta)-psi(sum(theta,4)));
        if sum(isnan(theta(:))) > 0
            flag=1; 
        end
        alpha = est_alpha(squeeze(sum(fM,1)),alpha);
        alpha(isnan(alpha)) = 0;
        alpha = alpha + .001;
        prob_diff = [prob_diff;[sum(abs(ao(:)-alpha(:)))/sum(ao(:)),mean(abs(to(:)-theta(:)))]];
        diff = sum(abs(ao(:)-alpha(:)))/sum(ao(:));
        if diff < .01 || iter >= max_iter
            conv = 1;
        end
    end

end

function C = const_mult(M,Z,A,pop)
    Z = Z+.001;
    C = zeros(size(M));
    for t=1:(size(Z,3)-1)
        for i=1:size(Z,1)
            for j=1:size(Z,2)
                s2 = .001;
                [N,th_n,th_p] = get_neighbors(i,j,A);
                for l=1:length(th_n)
                    s2 = s2 + M(N(l,1),N(l,2),t,th_p(l));
                end
                v = Z(i,j,t+1)-s2;

                for l=1:length(th_n)
                    C(N(l,1),N(l,2),t,th_p(l)) = v;
                end
            end
        end
    end

    C = 4*C/pop;
    C(C(:)>1) = 1;
    C(C(:)<-1) = -1;
    C = exp(C);
end

function alpha = est_alpha(M,alpha)   
    for i=1:1
        alpha_o = alpha;
        num = sum(psi(mean(mean(M+alpha,1),2))-psi(mean(mean(alpha,1),2)),3);
        den = sum(psi(mean(mean(sum(M+alpha,4),1),2))-psi(mean(mean(sum(alpha,4),1),2)),3);
        alpha = alpha.*(num./den);
        if mean(abs(alpha(:)-alpha_o(:))) < .01
            break;
        end
    end
end

function phi = init_phi(M,A)
    phi = zeros(size(M));
    for i=1:size(phi,1)
        for j=1:size(phi,2)
            [~,th_n,~] = get_neighbors(i,j,A);
            for t=1:size(phi,3)
                for n=1:length(th_n)
                    phi(i,j,t,th_n(n)) = 1;
                end
                %phi(i,j,t,5) = 100;
            end
        end
    end
end

function phi = normalize(phi,A)
    for t=1:(size(phi,3)-1)
        for i=1:size(phi,1)
            for j=1:size(phi,2)
                [~,th_n,~] = get_neighbors(i,j,A);
                val = sum(squeeze(phi(i,j,t,th_n)));
                phi(i,j,t,th_n) = phi(i,j,t,th_n)./val;
            end
        end
    end
end


function [mf,mb] = calc_messages(mf,mb,phi,A)
    %update messages foward
    for t=2:size(mf,3)
        for i=1:size(mf,1)
            for j=1:size(mf,2) 
                [N,th_n,th_p] = get_neighbors(i,j,A);
                mf_n = 0;
                for n=1:length(th_n)
                    mf_n = mf_n + phi(N(n,1),N(n,2),t-1,th_p(n))*mf(N(n,1),N(n,2),t-1);
                end
                mf(i,j,t) = mf_n;
                if isnan(mf_n)
                    flag=1;
                end
                %else
                %    mf_n = 1;
            end
        end
        if max(max(mf(:,:,t))) > 0
            mf(:,:,t) = mf(:,:,t)/max(max(mf(:,:,t)));
        else
            mf(:,:,t) = ones(size(mf,1),size(mf,2));
        end
    end
        
    for t=(size(mf,3)-1):-1:1
        for i=1:size(mf,1)
            for j=1:size(mf,2) 
                [N,th_n,~] = get_neighbors(i,j,A);
                mb_n=0;
                for n=1:length(th_n)
                    mb_n = mb_n + phi(i,j,t,th_n(n))*mb(N(n,1),N(n,2),t+1);
                end
                mb(i,j,t) = mb_n;
                if isnan(mb_n)
                    flag=1;
                end
            end
        end
        if max(max(mb(:,:,t))) > 0
            mb(:,:,t) = mb(:,:,t)/max(max(mb(:,:,t)));
        else
            mb(:,:,t) = ones(size(mb,1),size(mb,2));
        end
    end
end

function M = update_m(M,Z,mf,mb,phi,A)
    M = 0*M;
    for t=1:(size(M,3)-1)
        for i=1:size(M,1)
            for j=1:size(M,2)
                [N,th_n,~] = get_neighbors(i,j,A);
                for n=1:length(th_n)
                    if t > 1
                        M(i,j,t,th_n(n)) = phi(i,j,t,th_n(n))*mf(i,j,t)*mb(N(n,1),N(n,2),t+1);
                    else
                        M(i,j,t,th_n(n)) = phi(i,j,t,th_n(n))*mb(N(n,1),N(n,2),t+1);
                    end
                    if isnan(M(i,j,t,th_n(n)))
                        flag=1;
                    end
                end
                
                %Normalize
                s = sum(squeeze(M(i,j,t,:)));
                if s > 0
                    M(i,j,t,:) = Z(i,j,t)*M(i,j,t,:)/s;
                end
            end
        end
    end
end