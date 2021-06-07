function V = inv_alr(V)
    V = [V,zeros(size(V,1),1)];
    m = max(V,[],2);
    V = exp(V-m);
    V = V./sum(V,2);
end