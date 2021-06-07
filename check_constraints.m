function reg = check_constraints(Z,M,A)
    Zout = sum(M,4);
    Zin = zeros(size(Z));
    for t=1:(size(Z,3)-1)
        for i=1:size(Z,1)
            for j=1:size(Z,2)
                s2 = 0;
                [N,th_n,th_p] = get_neighbors(i,j,A);
                for l=1:length(th_n)
                    s2 = s2 + M(N(l,1),N(l,2),t,th_p(l));
                end
                Zin(i,j,t+1) = s2;
            end
        end
    end
    reg = abs(Zout(:,:,2:(end-1))-Zin(:,:,2:(end-1)));
    reg = mean(reg(:));
end