function save_temp_part2(folder_s2,T_s2,binT_s2,k_s2,X_s2,Y_s2,Z_s2,Vx_s2,Vy_s2,Vz_s2,Ax_s2,Ay_s2,Az_s2,Ntrack_s2)
    I=find(T_s2==binT_s2(k_s2));
    tp2.T=binT_s2(k_s2);
    %   X=X(I);
    %   Y=Y(I);
    tp2.X=X_s2(I);
    tp2.Y=Y_s2(I);
    tp2.Z=Z_s2(I);
    %   Theta=Theta(I);
    %   Rho=Rho(I);
    tp2.Vx=Vx_s2(I);
    tp2.Vy=Vy_s2(I);
    tp2.Vz=Vz_s2(I);
    tp2.Ax=Ax_s2(I);
    tp2.Ay=Ay_s2(I);
    tp2.Az=Az_s2(I);
    tp2.Ntrack=Ntrack_s2(I);

    save([folder_s2 filesep 'temp_part' filesep 'temp_part_' num2str(k) '.mat'],'tp2')
end