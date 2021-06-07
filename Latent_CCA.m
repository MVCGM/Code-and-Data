function [S,W1,W2,psi1,psi2,A,B,r] = Latent_CCA(phi1,phi2)

    [A,B,r] = canoncorr(phi1,phi2);
    
    mu1 = sum(phi1,1)/size(phi1,1);
    mu2 = sum(phi2,1)/size(phi2,1);
    
    P = diag(sqrt(r));
    I = eye(size(P,1));
    
    %Calculate CCA parameters and latent variable
    sig1 = cov(phi1);
    sig2 = cov(phi2);
    W1 = (sig1*A*P)';
    W2 = (sig2*B*P)';
    psi1 = sig1 - W1'*W1;
    psi2 = sig2 - W2'*W2;

    if abs(det(psi1)) < 10^-10 
        flag=1;
    end
    if abs(det(psi2)) < 10^-10 
        flag=1;
    end
    
    varS = I;
    muS = zeros(size(phi1,1),size(W1,1));
    
    V1 = W1*pinv(psi1);
    V2 = W2*pinv(psi2);
    
    S = (V1*W1'+V2*W2'+varS)\(V1*(phi1-mu1)'+V2*(phi2-mu2)'+varS*muS');
    S = S';
    
    if sum(sum(isnan(S)))>0
        flag=1;
    end
end