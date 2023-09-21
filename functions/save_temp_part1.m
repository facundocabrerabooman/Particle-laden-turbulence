function save_temp_part1(folder_s1,vtracks_s1,fields_s1,T_s1,Ntrack_s1,binT_s1,k_s1)
    tp1.T=binT_s1(k_s1);
    I_s1=find(T_s1==binT_s1(k_s1));
    tp1.Ntrack=Ntrack_s1(I_s1);
    
    for kfield_s1=1:numel(fields_s1)
        if isrow(vtracks_s1(1).(fields_s1{kfield_s1}))
            X_s1=[vtracks_s1.(fields_s1{kfield_s1})];
        else
            X_s1=vertcat(vtracks_s1.(fields_s1{kfield_s1}));
        end
        tp1.(fields_s1{kfield_s1})=X_s1(I_s1);
    end
   
    save([folder_s1 filesep 'temp_part' filesep 'temp_part_' num2str(k_s1) '.mat'],'tp1')
end