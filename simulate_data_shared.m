%Return true CGM models for both views
function [Z1,Z2,M1,M2,theta1,theta2,A] = simulate_data_shared(g1,g2,t,method)
    step=1;
    k=(2*step+1)^2-1;
    d=k;
    A = adj_matrix(g1,g2,step);

    p1 = k+1;
    p2 = p1;
    theta1n = zeros(g1,g2,t,p1);
    theta2n = zeros(g1,g2,t,p2); 
    
    id=1;
    for ts=1:t
        if mod(ts,g1)==1
            it1 = ceil(g1*rand(1));
            jt1 = ceil(g2*rand(1));
            it2 = ceil(g1*rand(1));
            jt2 = ceil(g2*rand(1));
        end
        for i=1:g1
            for j = 1:g2
                [N,idn] = get_neighbors(i,j,A);
                D1 = sqrt((it1-N(:,1)).^2+(jt1-N(:,2)).^2);
                D2 = sqrt((it2-N(:,1)).^2+(jt2-N(:,2)).^2);
                D = min(D1,D2);
                switch method
                    case 0
                        theta1n(i,j,ts,idn) = exp(-2*(D-d));
                        theta2n(i,j,ts,idn) = exp(-2*(D-d));
                    case 1
                        theta1n(i,j,ts,idn) = exp(-2*(D1-d));
                        theta2n(i,j,ts,idn) = exp(-2*(D2-d));
                    case 2
                        theta1n(i,j,ts,idn) = exp(-2*(D-d));
                        theta2n(i,j,ts,idn) = exp(2*(D-d));
                end
            end
            id = id+1;
        end
    end
    
    theta1 = zeros(g1,g2,t,p1);
    theta2 = zeros(g1,g2,t,p2); 
    for i=1:g1
        for j=1:g2
            [~,idn,~] = get_neighbors(i,j,A);
            theta1(i,j,:,idn) = theta1n(i,j,:,idn);
            theta2(i,j,:,idn) = theta2n(i,j,:,idn);
            theta1(i,j,:,:) = theta1(i,j,:,:)./sum(theta1(i,j,:,:),4);
            theta2(i,j,:,:) = theta2(i,j,:,:)./sum(theta2(i,j,:,:),4);
        end
    end
    
    Z1 = zeros(g1,g2,t);
    Z2 = Z1;
    M1 = zeros(g1,g2,t,p1);
    M2 = zeros(g1,g2,t,p2);
    
    %Seed population sources
    L1 = Z1(:,:,1);
    L2 = Z2(:,:,1);
    ids1 = datasample(1:length(L1(:)),2*g1,'Replace',false);
    ids2 = datasample(1:length(L2(:)),2*g1,'Replace',false);
    L1(ids1) = round(1000*rand(1,length(ids1)));
    L2(ids2) = round(1000*rand(1,length(ids2)));
    Z1(:,:,1) = L1;
    Z2(:,:,1) = L2;
    
    %Generate movement and population values for all times
    for k=1:(t-1)
        for i=1:g1
            for j=1:g2
                M1(i,j,k,:) = mnrnd(Z1(i,j,k),squeeze(theta1(i,j,k,:)));
                M2(i,j,k,:) = mnrnd(Z2(i,j,k),squeeze(theta2(i,j,k,:)));
                M1(i,j,k,:) = round(Z1(i,j,k)*squeeze(theta1(i,j,k,:)));
                M2(i,j,k,:) = round(Z2(i,j,k)*squeeze(theta2(i,j,k,:)));
                [N,idn,~] = get_neighbors(i,j,A);
                for p=1:length(idn)
                    Z1(N(p,1),N(p,2),k+1) = Z1(N(p,1),N(p,2),k+1)+M1(i,j,k,idn(p));
                    Z2(N(p,1),N(p,2),k+1) = Z2(N(p,1),N(p,2),k+1)+M2(i,j,k,idn(p));
                end
            end
        end
    end

end