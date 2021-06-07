function theta1 = matrix_form(theta,A)
    sz = size(theta);
    sz(3) = sz(3) -1;
    if sz(3) <= 0
        sz(3) = 1;
    end
    %convert into matrix forms
    theta1 = zeros(sz(1)*sz(2)*sz(3),sz(4));
    idx=1;
    for t=1:sz(3)
        for i=1:sz(1)
            for j=1:sz(2)
                [~,th_n,~] = get_neighbors(i,j,A);
                theta1(idx,th_n) = squeeze(theta(i,j,t,th_n))';
                idx = idx+1;
            end
        end
    end
end