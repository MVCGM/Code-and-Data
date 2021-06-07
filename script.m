function [mets,M] = script(reps,n,t,iters,method)
    M = [];
    mets = 0;
    ll = zeros(reps,3);
    Zd = ll;
    Td = ll;
    time = 0;
    for i = 1:reps
        i
        [m] = CGM_test(n,t,iters,method);
        M(:,i) = m;
        mets = mets+m(:,end);
    end
    mets = mets/reps;
    %save('M.mat','M');
end