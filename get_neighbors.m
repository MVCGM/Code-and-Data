%th_p are theta indicies from the previous time
%th_n are theta indicies going to the next time
function [N,idn,idp] = get_neighbors(i,j,A)
    g1 = size(A,1);
    g2 = size(A,2);
    N=[];
    idn = [];
    idp = [];
    
    len = max(max(squeeze(sum(sum(A,4),3))));
    step = (sqrt(len)-1)/2;
    
    id=0;
    for i2 = (i-step):(i+step)
        for j2 = (j-step):(j+step)
            id=id+1;
            if i2 <= 0 || i2 > g1
                continue;
            end
            if j2 <= 0 || j2 > g2
                continue;
            end
            N = [N;[i2,j2]];
            idn = [idn,id];
            idp = [idp,len-id+1];
        end
    end
end