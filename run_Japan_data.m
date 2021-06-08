function n = run_Japan_data(pref,Fout,max_iter)
    dir = "JapanPeopleFlow/"+pref+"/";
    switch pref
        case "Tokyo"
            days = ["07_01_10x10.mat","07_07_10x10.mat","10_07_10x10.mat","10_13_10x10.mat","12_16_10x10.mat","12_22_10x10.mat"];
        case "Osaka"
            days = ["08_08_10x10.mat","08_11_10x10.mat","09_16_10x10.mat","09_22_10x10.mat","12_24_10x10.mat","12_29_10x10.mat"];
        case "Nagoya"
            days = ["07_22_10x10.mat","07_28_10x10.mat","09_16_10x10.mat","09_22_10x10.mat","12_24_10x10.mat","12_29_10x10.mat"];
    end
    
    %Load data from file
    fM1 = [];
    fM2 = [];
    fZ1 = [];
    fZ2 = [];
    fpop1=[];
    fpop2=[];
    N = length(days);
    step1 = 1;%min(step_size(M1),1);
    step2 = 1;%min(step_size(M2),1);
    A1 = [];
    A2 = [];
    
    for i=1:N
        Fin = dir+days(i);
        L = load(Fin);
        M1 = L.M1;
        M2 = L.M2;
        if length(A1) == 0
            A1 = adj_matrix(size(M1,1),size(M1,2),step1);
            A2 = adj_matrix(size(M2,1),size(M2,2),step2);
        end
        M1 = convert_M(M1,step1,A1);
        M2 = convert_M(M2,step2,A2);
        fM1(i,:,:,:,:) = M1;
        fM2(i,:,:,:,:) = M2;
        fZ1(i,:,:,:) = sum(M1,4);
        fZ2(i,:,:,:) = sum(M2,4);
        fpop1(i) = sum(M1(:))/size(M1,3);
        fpop2(i) = sum(M2(:))/size(M2,3);
    end
    
    fM = fM1+fM2;
    tic
    [M,Z,phi,L] = Japan_CGM(ones(size(fM1)),max_iter,fpop1+fpop2,A1,fM1+fM2,fZ1+fZ2);
    t = toc;
    nae = 0;
    for i=1:size(M,1)
        nae = nae + Mdiff(squeeze(fM(i,:,:,:,:)),squeeze(M(i,:,:,:,:)));
    end
    nae = nae(1)/size(fM,1);
    save("base_all_"+Fout,"M","Z","t","L","nae");
    
    %Run algorithms, track time
    tic
    [M1,phi1,Z1,M2,phi2,Z2,L] = Japan_MVCGM(ones(size(fM1)),ones(size(fM2)),max_iter,fM1,fM2,fpop1,fpop2,A1,A2,fZ1,fZ2);
    t = toc;
    nae1 = 0;
    nae2 = 0;
    for i=1:size(M1,1)
        nae1 = nae1 + Mdiff(squeeze(fM1(i,:,:,:,:)),squeeze(M1(i,:,:,:,:)));
        nae2 = nae2 + Mdiff(squeeze(fM2(i,:,:,:,:)),squeeze(M2(i,:,:,:,:)));
    end
    nae1 = nae1(1)/size(M1,1);
    nae2 = nae2(1)/size(M2,1);
    save("multi_"+Fout,"M1","M2","Z1","Z2","t","L","nae1","nae2");
    flag=1
    
    tic
    [M1,Z1,phi1b,L1] = Japan_CGM(ones(size(fM1)),max_iter,fpop1,A1,fM1,fZ1);
    flag=2
    [M2,Z2,phi2b,L2] = Japan_CGM(ones(size(fM2)),max_iter,fpop2,A2,fM2,fZ2);
    t = toc;
    nae1 = 0;
    nae2 = 0;
    for i=1:size(M1,1)
        nae1 = nae1 + Mdiff(squeeze(fM1(i,:,:,:,:)),squeeze(M1(i,:,:,:,:)));
        nae2 = nae2 + Mdiff(squeeze(fM2(i,:,:,:,:)),squeeze(M2(i,:,:,:,:)));
    end
    nae1 = nae1(1)/size(fM1,1);
    nae2 = nae2(1)/size(fM2,1);
    save("base_"+Fout,"M1","M2","Z1","Z2","t","L1","L2","nae1","nae2");
    flag=3
    
    
    flag=4
    
    
    
    flag=3
    
    n=1;
end

function Mn = convert_M(M,step,A)
    k =(2*step+1)^2;
    Mn = zeros(size(M,1),size(M,2),size(M,3),k);
    for i=1:size(M,1)
        for j=1:size(M,2)
            for t=1:size(M,3)
                [N,th_n,~] = get_neighbors(i,j,A);
                for n=1:size(N,1)
                    Mn(i,j,t,th_n(n)) = M(i,j,t,N(n,1),N(n,2));
                end
            end
        end
    end
end