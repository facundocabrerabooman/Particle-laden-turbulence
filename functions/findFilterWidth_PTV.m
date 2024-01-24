function [s, m, w]=findFilterWidth_PTV(vel,field)

% [s, m, w]=findFilterWidth_PTV(vel,field);

w=1:170;
l=max(5,3*w);


for j=1:numel(w)
    %disp(sprintf('w = %i',w(j)));
    
    kerp = posfiltcoef(w(j),l(j));
    kerv = velfiltcoef(w(j),l(j));
    kera = accfiltcoef(w(j),l(j));
    
    xf=[];
    velxf=[];
    accxf=[];
    
    xf=cell2struct(arrayfun(@(X)(conv(X.(field),kerp,'valid')),vel,'UniformOutput',false),'xf');
    velxf=cell2struct(arrayfun(@(X)(conv(X.(field),kerv,'valid')),vel,'UniformOutput',false),'xf');
    accxf=cell2struct(arrayfun(@(X)(conv(X.(field),kera,'valid')),vel,'UniformOutput',false),'xf');
   
    L=arrayfun(@(x)(numel(x.xf)),xf);
    
    II = find(L~=0);
    
    m.x(j)=sum(arrayfun(@(x)(mean(x.xf)),xf(II)).*L(II))/sum(L(II));
    s.x(j)=sqrt(sum(arrayfun(@(x)(mean((x.xf-m.x(j)).^2)),xf(II)).*L(II))/sum(L(II)));
    m.vx(j)=sum(arrayfun(@(x)(mean(x.xf)),velxf(II)).*L(II))/sum(L(II));
    s.vx(j)=sqrt(sum(arrayfun(@(x)(mean((x.xf-m.vx(j)).^2)),velxf(II)).*L(II))/sum(L(II)));
    m.ax(j)=sum(arrayfun(@(x)(mean(x.xf)),accxf(II)).*L(II))/sum(L(II));
    s.ax(j)=sqrt(sum(arrayfun(@(x)(mean((x.xf-m.ax(j)).^2)),accxf(II)).*L(II))/sum(L(II)));
    
end

