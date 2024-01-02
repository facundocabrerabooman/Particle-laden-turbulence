% programme d'extraction des fragments de signal interessant% pour plusieurs voies d'acquisition simultanee%%% �2000 nmordant %%%%%%%%%%%%%% PARAMETRES%%%%%%%%%%%%%%%recherche des indicesseuil_haut=1e-3;seuil_bas=2e-5;seuil_long=300;% decodelong_supp=10000;dec_sup=60;% affichagept_sup=500;affiche_flag=1;%rate=10;%%%%%%%%%%%%%%%%%%%%%%%%%% RECHERCHE DES INDICES %%%%%%%%%%%%%%%%%%%%%%%%%% disp('Recherche des indices...');% preparation du signal...[bb2,aa2]=butter(6,2*100/1.28/span*rate);% if length(data1)==65536%     data1=data1(:,200:64500);% end;[NbVoies,long_tot]=size(data1);% recherche des maxima d'energietest=mean(abs(data1));q=find(test>seuil_haut);test(q)=seuil_haut;				% seuillage de l'energie avant filtragetest=filtfilt(bb2,aa2,test);	% filtrage de la fonction de test%plot(test);pauseq=find(test(400:(length(test)-200))>seuil_bas)+399;         % recherche des intervalles au dessus du seuilif ~isempty(q)    dq=diff(q);    q1=find(dq>1);    ind=zeros(length(q1)+1,2);    ind(1,1)=q(1);    nbfrag=1;    for i=1:length(q1)				% on raccorde les r�gions sous le seuil qui sont plus proches que seuil_long        if dq(q1(i))>seuil_long            ind(nbfrag,2)=q(q1(i));            ind(nbfrag+1,1)=q(q1(i)+1);            nbfrag=nbfrag+1;        end;    end;    ind(nbfrag,2)=q(length(q));    ind=ind(1:nbfrag,:);    disp(['Nombre total de fragments ',int2str(nbfrag)]);   else     nbfrag=0;    disp('pas de fragments!');end;clear q dq;%%%%%%%%%%%%% DECODAGE %%%%%%%%%%%%%disp(' ');disp('Extraction des vitesses...');%deb_ind=input('on commence par quel segment (n�) ?  [1]');%if isempty(deb_ind)%    deb_ind=1;%end;deb_ind=1;ind_fin_seg=0;nb_seg_fin=0;ind_deb_fin_seg=zeros(2,nbfrag);for z=deb_ind:nbfrag    disp(['segment n ',int2str(z),'/',int2str(nbfrag)]);    [mm,maxtest]=max(test(ind(z,1):ind(z,2)));		% max d'energie dans le segment=point de d�part    if (ind(z,1)+maxtest-1>ind_fin_seg)                % on commence par aller dans le sens du temps en partant du max        flag_init=1;        deb=(ind(z,1)+maxtest-1);        fin=min(long_tot,(ind(z,2)+long_supp));        if deb<dec_sup            deb=dec_sup+1;        end;        sig=data1(:,deb-dec_sup:fin);        MVA13mult1;        freq1=freq(1,:);        iHes1=iHes(1,:);        ind_fin_seg=deb-dec_sup+length(freq1)-1;        Puis1=Puis(1,:);                % on remonte le temps maintenant        sig=conj(data1(:,min(length(data1),deb+dec_sup):-1:max(1,(ind(z,1)-long_supp))));        flag_init=0;        NbS_0=NbStab(dec_sup);        Teta_0=freq(1:NbS_0,dec_sup)'*2*pi;        iHes_0=diag(iHes(1:NbS_0,dec_sup));        MVA13mult1;        freq=fliplr(freq);        freq=freq(1,:);        iHes=fliplr(iHes);        iHes=iHes(1,:);        ind_deb_seg=deb+dec_sup-length(freq)+1;        Puis=fliplr(Puis);        Puis=Puis(1,:);                % on construit une solution pour toute la fen�tre        if length(freq1)<2*dec_sup             freq1=[freq1,zeros(1,2*dec_sup)];            iHes1=[iHes1,zeros(1,2*dec_sup)];            Puis1=[Puis1,zeros(1,2*dec_sup)];        end;                freqres=[freq(1:length(freq)-3/2*dec_sup-1),0.5*freq(length(freq)-3/2*dec_sup:length(freq)-1/2*dec_sup)+0.5*freq1(dec_sup/2+1:3/2*dec_sup+1),freq1(3/2*dec_sup+2:length(freq1))];        iHesres=[iHes(1:length(iHes)-3/2*dec_sup-1),0.5*iHes(length(freq)-3/2*dec_sup:length(freq)-1/2*dec_sup)+0.5*iHes1(dec_sup/2+1:3/2*dec_sup+1),iHes1(3/2*dec_sup+2:length(freq1))];           Puisres=[Puis(1:length(Puis)-3/2*dec_sup-1),0.5*Puis(length(freq)-3/2*dec_sup:length(freq)-1/2*dec_sup)+0.5*Puis1(dec_sup/2+1:3/2*dec_sup+1),Puis1(3/2*dec_sup+2:length(freq1))];        Ltot=length(freq)+length(freq1)-dec_sup*2-1;        tres=linspace(0,(Ltot-1)/1.28/span*rate,Ltot);        if ind_fin_seg-ind_deb_seg>800            disp('		on le garde...')            nb_seg_fin=nb_seg_fin+1;            ind_deb_fin_seg(:,nb_seg_fin)=[ind_deb_seg;ind_fin_seg];            com=['seg',int2str(nb_seg_fin),'=[freqres;iHesres;Puisres];'];            eval(com);        else            disp('		trop court');        end;        ind_deb_fin_seg=ind_deb_fin_seg(:,1:nb_seg_fin);                % affichages divers        if affiche_flag==1            [B,F,T]=specgram(data1(1,max(1,ind_deb_seg-pt_sup):min(long_tot,ind_fin_seg+pt_sup)),256,1.28*span/rate,64,57);B=imshift(B);            dec_tempB=32/1.28/span*rate;            dec_tempsup=(dec_sup+1)/1.28/span*rate;            Dt=pt_sup/1.28/span*rate;            subplot(311);            imagesc(T+dec_tempB,F-F(129),abs(B).^0.1);            hold on;            plot(length(freq)/1.28/span*rate-dec_tempsup+Dt,0,'k*');            plot(tres+Dt,(freqres)'*1.28*span/rate,'b');            title([fname,'   rack 1   segment n�',int2str(z),'   indice de debut ',int2str(ind(z,1)+maxtest-1)]);            hold off;            subplot(312);            hold off;            plot(tres+Dt,sqrt(iHesres(1,:)'));            axis([dec_tempB,max(T+dec_tempB),0,0.1]);hold off;            subplot(313);            hold off;            plot(tres+Dt,sqrt(Puisres(1,:)));            axis([dec_tempB,max(T+dec_tempB),0,max(sqrt(Puisres(1,:)))]);hold off;            pause(0.1);        end;    end;   end;clear aa2 bb2 N O a i dec_sup long_supp pt_sup seuil_haut test long_tot q1 seuil_long seuil_basclear B Teta_0 freq1 maxtest Dt mm F com Ltot iHes1 iHes iHes_0 trG NbS_0 debclear NbStab deb_ind ind Puis dec_tempB Puis1 dec_tempsup ind_deb_seg Puisres tresclear T freqres nbfrag iHesres ind_fin_seg freq fin z sig NbVoiesclear affiche_flag b f ind_Puis n t