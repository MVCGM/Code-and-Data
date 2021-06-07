function [M,phi,Ze,LM] = Vec_CGM(M,max_iter,tM,pop,A,Z,bias)
    %Construct initial phi and messages
    phi = init_phi(M,A);
    alpha = 2*ones(size(M,1),size(M,2));
    
    mf = ones(size(M,1),size(M,3));
    mb = mf;

    theta = rand(size(phi));
    theta = theta./sum(theta,2);

    conv = 0;

    %Loop till convergence 
    iter = 0;
    LM = zeros(max_iter,2);
    beta = 0.9;

    while conv == 0
        iter = iter+1;
        
        C = const_mult(M,Z,pop);
        
        phi = theta.*bias.*C;
        
        %Calculate messages
        [mf,mb] = calc_messages(mf,mb,phi,A);
        
        %Update M
        M_n = update_m(M,Z,mf,mb,phi,A);
        
        M_n = beta*M + (1-beta)*M_n;
        
        diff = sum(abs(M_n(:)-M(:)))/sum(M(:));

        if isnan(diff)
            flag=1;
            break;
        end
        
        M = M_n;
        
        acc = Mdiff(tM,M_n);

        LM(iter,:) = [acc(1),diff];        
            
        if mean(diff) < .005 || iter >= max_iter
            conv=1;
        end
        
        theta = alpha+M;
        theta = exp(psi(theta)-psi(sum(theta,2)));
        alpha = est_alpha(M,alpha);

    end
    Ze = mf.*mb;
    for t=1:size(M,3)
        if sum(sum(Ze(:,t))) == 0
            Ze(:,t) = 1;
        end
        Ze(:,t) = Ze(:,t) * pop(t) /sum(Ze(:,t));
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

function alpha = est_alpha(M,alpha)   
    for i=1:1
        alpha_o = alpha;
        num = sum(psi(M+alpha)-psi(alpha),3);
        den = sum(psi(mean(M+alpha,2))-psi(mean(alpha,2)),3);
        alpha = alpha.*(num./den);
        alpha(isnan(alpha)) = 0;
        alpha(alpha(:) < .0001) = .0001;
        if mean(abs(alpha(:)-alpha_o(:))) < .01
            break;
        end
    end
end


function phi = init_phi(M,A)
    phi = zeros(size(M));
    for t=1:size(phi,3)
        phi(:,:,t) = A;
    end
            
end


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