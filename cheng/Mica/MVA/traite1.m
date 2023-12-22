for iii=6:6
    for jjj=2:60
        %if ((iii==3)&jjj>20)|iii>3
            fname=['run',int2str(iii),'_',int2str(jjj)];
            disp(fname);
            if exist(fname,'file')
                traite2;
                %cmd=['load run',int2str(iii),'_',int2str(jjj)];
                %eval(cmd);
                %clear seg*
                %clear ind_deb_fin_seg
                extract13mult1;
                if nb_seg_fin>0
                    cmd=['save run',int2str(iii),'_',int2str(jjj),'_1 data1 nb_seg_fin ind_deb_fin_seg seg* span rate'];
                    eval(cmd);
                    clear seg* ind_deb_fin_seg nb_seg_fin;
                end;
                
                extract13mult2;
                if nb_seg_fin>0  
                    cmd=['save run',int2str(iii),'_',int2str(jjj),'_2 data2 nb_seg_fin ind_deb_fin_seg seg* span rate'];
                    eval(cmd);
                    clear data* seg* ind_deb_fin_seg nb_seg_fin ;
                end;
            end;
            %end;
    end;
end;





