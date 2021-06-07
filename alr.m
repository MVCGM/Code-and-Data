function V = alr(V)
    min_v = .000001;
    V(V < min_v) = min_v;
    V = log(V(:,1:(end-1))./V(:,end));
end