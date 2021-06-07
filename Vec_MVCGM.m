function [M1,phi1,Ze1,M2,phi2,Ze2,LM] = Vec_MVCGM(M1,M2,max_iter,tM1,tM2,pop1,pop2,A1,A2,Z1,Z2,bias)
    %Construct initial phi and messages
    phi1 = init_phi(M1,A1);
    phi2 = init_phi(M2,A2);
    
    mf1 = ones(size(M1,1),size(M1,3));
    mf2 = ones(size(M2,1),size(M2,3));
    mb1 = mf1;
    mb2 = mf2;

    theta1 = rand(size(phi1));
    theta1 = theta1./sum(theta1,2);
    theta2 = rand(size(phi2));
    theta2 = theta2./sum(theta2,2);

    conv = 0;

    %Loop till convergence 
    iter = 0;
    LM = zeros(max_iter,4);
    beta = 0.9;

    while conv == 0
        iter = iter+1;
        C1 = const_mult(M1,Z1,pop1);
        C2 = const_mult(M2,Z2,pop2);
        
        phi1 = theta1.*bias.*C1;
        phi2 = theta2.*bias.*C2;
        
        %Calculate messages
        [mf1,mb1] = calc_messages(mf1,mb1,phi1,A1);
        [mf2,mb2] = calc_messages(mf2,mb2,phi2,A2); 
        
        %Update M
        M1_n = update_m(M1,Z1,mf1,mb1,phi1,A1);
        M2_n = update_m(M2,Z2,mf2,mb2,phi2,A2);
        
        M1_n = beta*M1 + (1-beta)*M1_n;
        M2_n = beta*M2 + (1-beta)*M2_n;
        
        diff = sum(abs(M1_n(:)-M1(:)))/sum(M1(:))+sum(abs(M2_n(:)-M2(:)))/sum(M2(:));

        if isnan(diff)
            flag=1;
        end
        
        M1 = M1_n;
        M2 = M2_n;
        
        acc1 = Mdiff(tM1,M1_n);
        acc2 = Mdiff(tM2,M2_n);

        LM(iter,:) = [acc1(1),acc2(1),acc1(1)+acc2(1),diff];        
            
        if mean(diff) < .005 || iter >= max_iter
            conv=1;
        end
        
        [theta1n,theta2n] = multiview_thetas(M1,M2);
        diff = .5*sum(abs(theta1n(:)-theta1(:)))/sum(theta1(:)) + .5*sum(abs(theta2n(:)-theta2(:)))/sum(theta2(:));
        if diff < .005 || iter >= max_iter
            %conv=1;
        end
        b = 0;
        theta1 = b*theta1+(1-b)*theta1n;
        theta2 = b*theta2+(1-b)*theta2n;

    end
    Ze1 = mf1.*mb1;
    Ze2 = mf2.*mb2;
    for t=1:size(M1,3)
        if sum(sum(Ze1(:,t))) == 0
            Ze1(:,t) = 1;
        end
        if sum(sum(Ze2(:,t))) == 0
            Ze2(:,t) = 1;
        end
        Ze1(:,t) = Ze1(:,t) * pop1(t) /sum(Ze1(:,t));
        Ze2(:,t) = Ze2(:,t) * pop2(t) /sum(Ze2(:,t));
    end
end

function C = const_mult(M,Z,pop)
    Z = Z+.001;
    C = zeros(size(M));
    for t=1:(size(Z,2)-1)
        for i=1:size(Z,1)
                s2 = .001;
                for j=1:size(M,2)
                    s2 = s2 + M(j,i,t);
                end
                v = Z(i,t+1)-s2;
                v=4*v/pop(t);
                for j=1:size(M,2)
                    C(j,i,t) = v;
                end
        end
    end
    C(C(:)>1) = 1;
    C(C(:)<-1) = -1;
    C = exp(C);
end

function [T1,T2] = multiview_thetas(fM1,fM2)
    T1 = zeros(size(fM1));
    T2 = zeros(size(fM2));

    M1 = sum(fM1,3);
    M2 = sum(fM2,3);
    
    %alr transform on M values
    phi1 = alr(M1);
    phi2 = alr(M2);

    %Estimate latent variables under normal form
    mu1 = mean(phi1);
    mu2 = mean(phi2);
        
    try
        [S,W1,W2,Psi1,Psi2,~,~,r] = Latent_CCA(phi1,phi2);
    catch
        flag=1;
    end
    
    %estimate sum(alpha) of Psi1 and Psi2
    p1 = cov2alpha(Psi1);
    p2 = cov2alpha(Psi2);
    
    %Compute expected values
    E1 = mu1+S*W1;
    E2 = mu2+S*W2;
    alpha1 = inv_alr(E1);
    alpha2 = inv_alr(E2);
    
    for t=1:size(fM1,3)
        for i=1:size(fM1,1)
            Tt = p1*alpha1(i,:)+fM1(i,:,t);
            Tt = exp(psi(Tt)-psi(sum(Tt)));
            Tt(isnan(Tt)) = 0;
            T1(i,:,t) = Tt;
        end
    end
    
    for t=1:size(fM2,3)
        for i=1:size(fM2,1)
            Tt = p2*alpha2(i,:)+fM2(i,:,t);
            Tt = exp(psi(Tt)-psi(sum(Tt)));
            Tt(isnan(Tt)) = 0;
            T2(i,:,t) = Tt;
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
    for t=1:size(phi,3)
        phi(:,:,t) = A;
    end
            
end

%May need to play with message passing order here
function [mf,mb] = calc_messages(mf,mb,phi,A)
    ms = size(mf);
    min_v = 10^-300;
    %Forward message pass
    for t=1:ms(2)
        if t>1
            for i=1:ms(1)
                mfn = 0;
                for j=1:ms(1)
                    if A(j,i) == 1
                        mfn = mfn + phi(j,i,t)*mf(j,t-1);
                    end
                end
                mf(i,t) = mfn;
            end
        end
        mf(:,t) = message_norm(mf(:,t),min_v);
    end

    %Backward message pass
    for t=(ms(2)):-1:1
        if t < ms(2)
            for i=1:ms(1)
                mbn = 0;
                for j=1:ms(1)
                    if A(i,j) == 1
                        mbn = mbn + phi(i,j,t)*mb(j,t+1);
                    end
                end
                mb(i,t) = mbn;
            end
        end
        mb(:,t) = message_norm(mb(:,t),min_v);
    end
end

function [mt] = message_norm(mt,mv)
    mm = max(mt(:));
    if mm > 0
        mt = mt/max(mt(:));
    else
        mt(:) = 1;
    end
    if min(mt(:)) == 0
        flag=1;
    end
    mt(mt(:) == 0) = mv;
end

function M = update_m(M,Z,mf,mb,phi,A)
    M = 0*M;
    for t=1:(size(M,3)-1)
        for i=1:size(M,1)
            for j=1:size(M,2)
                if A(i,j)==1
                    M(i,j,t) = phi(i,j,t)*mf(i,t)*mb(j,t+1);
                end
            end
            %Normalize
            s = sum(M(i,:,t));
            if s > 0
                M(i,:,t) = Z(i,t)*M(i,:,t)/s;
            end
        end
    end
end