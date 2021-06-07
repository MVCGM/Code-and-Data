function theta = inv_matrix_form(M,theta)
    sz = size(theta);
    sz(3) = sz(3) -1;
    if sz(3) <= 0
        sz(3) = 1;
    end
    
    idx=1;
    for t=1:sz(3)
        for i=1:sz(1)
            for j=1:sz(2)
                theta(i,j,t,:) = M(idx,:);
                idx = idx+1;
            end
        end
    end
end