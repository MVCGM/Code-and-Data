function mets = CGM_test(n,t,max_iter,method)
    %Simulate data with two views
    [Z1,Z2,M1,M2,theta1,theta2,A] = simulate_data_shared(n,n,t,method);

    pop1 = sum(Z1(:))/size(Z1,3);
    pop2 = sum(Z2(:))/size(Z2,3);
    

    Mx = zeros(size(M1));

    fpop1(1) = pop1;
    fpop2(1) = pop2;
    fMx(1,:,:,:,:) = Mx;
    fM1(1,:,:,:,:) = M1;
    fM2(1,:,:,:,:) = M2;
    fZ1(1,:,:,:) = Z1;
    fZ2(1,:,:,:) = Z2;

    time = [0;0];
    
    %Run algorithms
    tic;
    [Mp1,~,Zp1,Mp2,~,Zp2,Lp] = Sim_MVCGM(fMx,fMx,max_iter,fM1,fM2,fpop1,fpop2,A,A,fZ1,fZ2);
    time(2) = toc;
    Lp = squeeze(Lp);
    Mp1 = squeeze(Mp1);
    Mp2 = squeeze(Mp2);
    Zp1 = squeeze(Zp1);
    Zp2 = squeeze(Zp2);
    MDp1 = Mdiff(M1,Mp1);
    MDp2 = Mdiff(M2,Mp2);
    MDpc = Mdiff(M1+M2,Mp1+Mp2);
    Zdist7 = Mdiff(Z1,Zp1);
    Zdist7 = Zdist7(1);
    Zdist8 = Mdiff(Z2,Zp2);
    Zdist8 = Zdist8(1);
    cnp = check_constraints(Zp1,Mp1,A)+check_constraints(Zp2,Mp2,A);
    
    
    tic;
    [Ms1,Zs1,~,Ls1] = Sim_CGM(fMx,max_iter,fpop1,A,fM1,fZ1);
    [Ms2,Zs2,~,Ls2] = Sim_CGM(fMx,max_iter,fpop2,A,fM2,fZ2);
    time(1) = toc;
    Ls1 = squeeze(Ls1);
    Ls2 = squeeze(Ls2);
    Ms1 = squeeze(Ms1);
    Ms2 = squeeze(Ms2);
    Zs1 = squeeze(Zs1);
    Zs2 = squeeze(Zs2);
    MDs1 = Mdiff(M1,Ms1);
    MDs2 = Mdiff(M2,Ms2);
    MDsc = Mdiff(M1+M2,Ms1+Ms2);
    Zdist1 = Mdiff(Z1,Zs1);
    Zdist1 = Zdist1(1);
    Zdist2 = Mdiff(Z2,Zs2);
    Zdist2 = Zdist2(1);
    cns = check_constraints(Zs1,Ms1,A)+check_constraints(Zs2,Ms2,A);
    
    

    mets = [cns;cnp;MDs1(1);MDs2(1);MDsc(1);MDp1(1);MDp2(1);MDpc(1);time];
end

