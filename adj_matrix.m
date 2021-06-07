function A = adj_matrix(g1,g2,step)
    %defines adjacencies by number of 'steps' away they can go.
    %diagonal is considered one step
    A = zeros(g1,g2,g1,g2);
    
    for i=1:g1
        for j=1:g2
            for i2 = (i-step):(i+step)
                if i2 <= 0 || i2 > g1
                    continue;
                end
                for j2 = (j-step):(j+step)
                    if j2 <= 0 || j2 > g2
                        continue;
                    end
                    A(i,j,i2,j2) = 1;
                end
            end
        end
    end
    
end