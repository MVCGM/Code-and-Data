function n = run_US_data(Fout,max_iter)
    dir = "US Flow/";
    
    %Load data from file
    L = load(dir+"JobM.mat");
    M1 = L.M;
    L = load(dir+"MobM.mat");
    M2 = L.M;
    
    T = readtable(dir+"Distances.csv");
    P = table2array(T(:,[3,4]));
    D = zeros(size(P,1),size(P,1));
    
    for i=1:size(D,1)
        for j=1:size(D,2)
            if i==j
                v = 0;
            else
                v = circle_dist(P(i,:),P(j,:));
            end
            D(i,j) = v;
        end
    end
    sig = 100000;
    D = exp(-D.^2/sig);
    
    %Remove unsupported states
    ids = ~isnan(squeeze(sum(M1(:,1,:),3)));
    M1 = M1(ids,ids,:);
    M2 = M2(ids,ids,:);
    D = D(ids,ids);
    
    A1 = ones(size(M1,1),size(M1,2));
    A2 = ones(size(M2,1),size(M2,2));

    tZ1 = squeeze(sum(M1,2));
    tZ2 = squeeze(sum(M2,2));
    
    pop1 = sum(tZ1,1);
    pop2 = sum(tZ2,1);
    
    tM1 = M1;
    tM2 = M2;
    
    tic
    [M1,phi1,Z1,M2,phi2,Z2,L] = Vec_MVCGM(ones(size(M1)),ones(size(M2)),max_iter,tM1,tM2,pop1,pop2,A1,A2,tZ1,tZ2,D);
    t = toc;
    nae1 = Mdiff(tM1,M1);
    nae1 = nae1(1);
    nae2 = Mdiff(tM2,M2);
    nae2 = nae2(1);
    save("multi_"+Fout,"M1","M2","Z1","Z2","t","L","nae1","nae2");
    flag=1
    
    tic
    [M1,Z1,phi1,L1] = Vec_CGM(ones(size(M1)),max_iter,tM1,pop1,A1,tZ1,D);
    flag=2
    [M2,Z2,phi2,L2] = Vec_CGM(ones(size(M1)),max_iter,tM2,pop2,A2,tZ2,D);
    t = toc;
    nae1 = Mdiff(tM1,M1);
    nae1 = nae1(1);
    nae2 = Mdiff(tM2,M2);
    nae2 = nae2(1);
    save("base_"+Fout,"M1","M2","Z1","Z2","t","L1","L2","nae1","nae2");
    
    flag=3
    
    
    n=1;
end

function d = circle_dist(p1,p2)
    p1 = deg2rad(p1);
    p2 = deg2rad(p2);
    lat1 = p1(1);
    lat2 = p2(1);
    long1 = p1(2);
    long2 = p2(2);
    
    d = 3963.0 * acos((sin(lat1)*sin(lat2)) + cos(lat1)*cos(lat2)*cos(long2 - long1));
end