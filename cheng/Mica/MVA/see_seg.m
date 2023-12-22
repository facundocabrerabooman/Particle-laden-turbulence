[ff]=see_seg(vel,ideb,iend,p);

%(figure)=see_seg(vel,ideb,iend,pause);
figure; hold on;
for i=ideb:iend;
    plot(abs(vel.data(vel.good(i)).seg).^2);
    text(150,7e-7,num2str(i))
    pause(p);
    clear text only
end



