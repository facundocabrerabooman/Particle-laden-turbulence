function [Lt,La,alpha,xcentre,ycentre,Rcp]=volumeacoustique(angle,frequence,d,diametre)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
% SYNTAXE : [Lt,La]=diffraction(angle,frequence,distance,diametre) %
%                                                                  %
%                                                                  %
% Ce script permet de calculer les deux longueurs caractéristiques %
% Lt et La du volume de mesure en tenant compte de la diffraction  %
% des capteurs. Pour ce calcul, il est nécessaire de connaître les %
% paramètres suivants :                                            %
% 	angle : angle de diffusion (par rapport au vis à vis)      %
% 	frequence : fréquence d'émission                           %
% 	d : distance entre les deux capteurs                        %
%       diametre : diametre des capteurs                           %
%   pa=distance : demi-parcours acoustique                         %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Calcul de l'angle de diffraction dteta
c=340;
dteta=asin(3.83*340./(pi*frequence*diametre))
%dteta=0;
teta=angle*pi/180
pa=d/(2*cos(teta))
Rcp=((diametre.^2*frequence/4/c)-(c/4/frequence))
if 2*pa*tan(dteta)<diametre
    La=diametre./cos(teta);
    Lt=La./tan(teta);
    xcentre=La/2;
    ycentre=La/2;
else
    La=2*pa.*cos(teta).*tan(dteta).*(cos(teta)+sin(teta).*tan(teta+dteta))./cos(teta-dteta);
    Lt=2*La./(tan(teta+dteta)+tan(teta-dteta));
    xcentre=0.5.*Lt.*tan(teta+dteta);
    ycentre=0.5.*Lt.*tan(teta-dteta);
end

alpha=dteta*180/pi;