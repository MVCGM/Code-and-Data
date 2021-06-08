function [fM1,fphi1,fZe1,fM2,fphi2,fZe2,LM] = Japan_MVCGM(fM1,fM2,max_iter,ftM1,ftM2,fpop1,fpop2,A1,A2,fZ1,fZ2)
    N = size(fM1,1);

    fphi1 = ones(size(fM1));
    fphi2 = ones(size(fM2));

    fmf1 = ones(N,size(fM1,2),size(fM1,3),size(fM1,4));
    fmf2 = ones(N,size(fM2,2),size(fM2,3),size(fM2,4));
    fmb1 = fmf1;
    fmb2 = fmf2;
    parfor i=1:N
        M1 = squeeze(fM1(i,:,:,:,:));
        M2 = squeeze(fM2(i,:,:,:,:));
        Z1 = squeeze(fZ1(i,:,:,:));
        Z2 = squeeze(fZ2(i,:,:,:));
        mf1 = squeeze(fmf1(i,:,:,:));
        mb1 = squeeze(fmb1(i,:,:,:));
        mf2 = squeeze(fmf2(i,:,:,:));
        mb2 = squeeze(fmb2(i,:,:,:));
        phi1 = init_phi(M1,A1);
        phi2 = init_phi(M2,A2);
        M1 = update_m(M1,Z1,mf1,mb1,phi1,A1);
        M2 = update_m(M2,Z2,mf2,mb2,phi2,A2);
        C1 = const_mult(M1,Z1,A1,fpop1(i));
        C2 = const_mult(M2,Z2,A2,fpop2(i));
        phi1 = phi1.*C1;
        phi2 = phi2.*C2;
        fM1(i,:,:,:,:) = update_m(M1,Z1,mf1,mb1,squeeze(fphi1(i,:,:,:,:)),A1);
        fM2(i,:,:,:,:) = update_m(M2,Z2,mf2,mb2,squeeze(fphi2(i,:,:,:,:)),A2);
        fphi1(i,:,:,:,:) = phi1;
        fphi2(i,:,:,:,:) = phi2;
    end
    theta1 = normalize(rand(size(squeeze(fphi1(1,:,:,:,:)))),A1);
    theta2 = normalize(rand(size(squeeze(fphi2(1,:,:,:,:)))),A2);
    
    conv = 0;
    
    %Loop till convergence 
    iter = 0;
    LM = zeros(N,max_iter,8);
    prob_diff=[];
    beta = 0.5;
    while conv == 0
        iter = iter+1;
        %For each dataset...
        diff = ones(1,N);

        for f=1:N
            mf1 = squeeze(fmf1(f,:,:,:));
            mb1 = squeeze(fmb1(f,:,:,:));
            mf2 = squeeze(fmf2(f,:,:,:));
            mb2 = squeeze(fmb2(f,:,:,:));
            M1 = squeeze(fM1(f,:,:,:,:));
            M2 = squeeze(fM2(f,:,:,:,:));
            Z1 = squeeze(fZ1(f,:,:,:));
            Z2 = squeeze(fZ2(f,:,:,:));
            pop1 = squeeze(fpop1(f));
            pop2 = squeeze(fpop2(f));
            tM1 = squeeze(ftM1(f,:,:,:,:));
            tM2 = squeeze(ftM2(f,:,:,:,:));
        
            C1 = const_mult(M1,Z1,A1,pop1);
            C2 = const_mult(M2,Z2,A2,pop2);
            phi1 = theta1.*C1;
            phi2 = theta2.*C2;
        
            %Calculate messages
            [mf1,mb1] = calc_messages(mf1,mb1,phi1,A1);
            [mf2,mb2] = calc_messages(mf2,mb2,phi2,A2);
            

        
            Ze1 = mf1.*mb1;
            Ze2 = mf2.*mb2;
            for t=1:size(M1,3)
                if sum(sum(Ze1(:,:,t))) == 0
                    Ze1(:,:,t) = 1;
                end
                if sum(sum(Ze2(:,:,t))) == 0
                    Ze2(:,:,t) = 1;
                end
                Ze1(:,:,t) = Ze1(:,:,t) * pop1 /sum(sum(Ze1(:,:,t)));
                Ze2(:,:,t) = Ze2(:,:,t) * pop2 /sum(sum(Ze2(:,:,t)));
            end
        
            %Update M
            M1_n = update_m(M1,Z1,mf1,mb1,phi1,A1);
            M2_n = update_m(M2,Z2,mf2,mb2,phi2,A2);
        
        
            M1_n = beta*M1 + (1-beta)*M1_n;
            M2_n = beta*M2 + (1-beta)*M2_n;
        
            diff(f) = sum(abs(M1_n(:)-M1(:)))/sum(M1(:))+sum(abs(M2_n(:)-M2(:)))/sum(M2(:));
        
            if isnan(diff(f))
                flag=1;
            end
        
            M1 = M1_n;
            M2 = M2_n;

        
            acc1 = Mdiff(tM1,M1_n);
            acc2 = Mdiff(tM2,M2_n);
            acc(f,:) = [acc1(1),acc2(1)];
        
            c1 = check_constraints(Ze1,M1,A1);
            c2 = check_constraints(Ze2,M2,A2);
            ct = sum(c1+c2);

            LM(f,iter,:) = [c1,c2,acc1(1),acc2(1),acc1(1)+acc2(1),diff(f),sum(sum(sum(abs(Z1-Ze1)))),sum(sum(sum(abs(Z2-Ze2))))];

        

            fmf1(f,:,:,:) = mf1;
            fmb1(f,:,:,:) = mb1;
            fmf2(f,:,:,:) = mf2;
            fmb2(f,:,:,:) = mb2;
            fM1(f,:,:,:,:) = M1;
            fM2(f,:,:,:,:) = M2;
            fZe1(f,:,:,:) = Ze1;
            fZe2(f,:,:,:) = Ze2;
            
        end

        [theta1n,theta2n] = est_theta(fM1,fM2);
        %diff = .5*sum(abs(theta1n(:)-theta1(:)))/sum(theta1(:)) + .5*sum(abs(theta2n(:)-theta2(:)))/sum(theta2(:));
        diff = mean(diff);
        if diff < .0025 || iter >= max_iter
            conv=1;
        end

        %LM(1,iter,6) = diff;
        theta1 = theta1n;
        theta2 = theta2n;
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

function M = matrix_convert(fM)
    M = squeeze(sum(sum(fM,1),2));
end

function [T1,T2,W1,W2] = est_theta(fM1,fM2)
    fM1 = squeeze(sum(fM1,1));
    fM2 = squeeze(sum(fM2,1));
    T1 = zeros(size(fM1));
    T2 = zeros(size(fM2));
    %transform all M values into a 2D matrix form
        M1 = matrix_convert(fM1(:,:,:,:));
        M2 = matrix_convert(fM2(:,:,:,:));
    
        n = size(M1,1);
        %alr transform on M values
        phi1 = alr(M1);
        phi2 = alr(M2);
    
        %Estimate latent variables under normal form
        mu1 = mean(phi1);
        mu2 = mean(phi2);
        
    
        try
        [S,W1,W2,Psi1,Psi2,U1,U2,r] = Latent_CCA(phi1,phi2);
        r(r > .999) = .999;
        catch
            flag=1;
        end
        
        if sum(isnan(S(:))) > 0
            flag=1;
        end
    
        %Compute expected values
        E1 = S*W1+mu1;
        E2 = S*W2+mu2;
        alpha1 = inv_alr(E1);
        alpha2 = inv_alr(E2);
        
        %estimate sum(alpha) of Psi1 and Psi2
        p1 = cov2alpha(Psi1);
        p2 = cov2alpha(Psi2);

        M1 = fM1;
        M2 = fM2;
        id=1;

        
        for i=1:size(M1,1)
            for j=1:size(M1,2)
                for t=1:size(M1,3)
                    Tt = p1*alpha1(t,:)+squeeze(M1(i,j,t,:))';
                    Tt = exp(psi(Tt)-psi(sum(Tt)));
                    Tt(isnan(Tt)) = 0;
                    T1(i,j,t,:) = Tt;
                    id=id+1;
                end
                
            end
        end
        id=1;
        
        for i=1:size(M2,1)
            for j=1:size(M2,2)
                for t=1:size(M2,3)
                    Tt = p2*alpha2(t,:)+squeeze(M2(i,j,t,:))';
                    Tt = exp(psi(Tt)-psi(sum(Tt)));
                    Tt(isnan(Tt)) = 0;
                    T2(i,j,t,:) = Tt;
                    id=id+1;
                end
                
            end
        end
end

function sp = cov2alpha(C)
    k = size(C,1)+1;
    ca = (sum(C(:))-sum(diag(C)))/(size(C,1)^2-size(C,1));
    if ca < 0 || ca > min(diag(C))
        ca = 0;
    end
    palpha = zeros(1,k);
    if ca > 0
        palpha(k) = inv_trigam(ca);
    end
    for i=1:k-1
        if (C(i,i)-ca) > 0
            palpha(i) = inv_trigam(C(i,i)-ca);
        end
    end
    sp = sum(palpha);
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
                if isnan(val)
                    flag=1;
                end
                if val > 0
                    phi(i,j,t,th_n) = phi(i,j,t,th_n)./val;
                end
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