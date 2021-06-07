function d = Mdiff(M1,M2)
    d = 0;
    s = 0;
    p = 0;
    ms = 0;
    md = 0;
    mp = 0;
    for t=1:(size(M1,3)-1)
        for i=1:size(M1,1)
            for j=1:size(M1,2)
                th_n = 1:size(M1,4);
                for n=1:length(th_n)
                    d = d + abs(M1(i,j,t,th_n(n)) - M2(i,j,t,th_n(n)));
                    s = s + M1(i,j,t,th_n(n));
                end
                [mv1,mid1] = max(M1(i,j,t,th_n));
                [mv2,mid2] = max(M2(i,j,t,th_n));
                if mid1 == mid2
                    p = p + 1;
                    md = md + abs(mv1-mv2);
                else
                    md = md + abs(mv2 - M1(i,j,t,th_n(mid2)));
                end
                mp = mp+1;
                ms = ms + M1(i,j,t,th_n(mid2));
            end
        end
    end
    
    d = d /s;
    md = md/ms;
    p = p/mp;
    d = [d;md;p];
end