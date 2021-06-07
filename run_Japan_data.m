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
    
    %Run algorithms, track time
    tic
    [M1m,phi1m,Z1m,M2m,phi2m,Z2m,Lm] = Japan_MVCGM(ones(size(fM1)),ones(size(fM2)),max_iter,fM1,fM2,fpop1,fpop2,A1,A2,fZ1,fZ2);
    tm = toc;
    save("multi_"+Fout,"M1m","M2m","Z1m","Z2m","tm","Lm");
    flag=1
    
    tic
    [M1b,Z1b,phi1b,Lb1] = Japan_CGM(ones(size(fM1)),max_iter,fpop1,A1,fM1,fZ1);
    flag=2
    [M2b,Z2b,phi2b,Lb2] = Japan_CGM(ones(size(fM2)),max_iter,fpop2,A2,fM2,fZ2);
    tb = toc;
    save("base_"+Fout,"M1b","M2b","Z1b","Z2b","tb","Lb1","Lb2");
    flag=3
    
    tic
    [Mb,Zb,phib,Lb] = Japan_CGM(ones(size(fM1)),max_iter,fpop1+fpop2,A1,fM1+fM2,fZ1+fZ2);
    tt = toc;
    save("base_all_"+Fout,"Mb","Zb","tt","Lb");
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