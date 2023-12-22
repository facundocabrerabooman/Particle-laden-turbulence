function MSD = LagrangianStatMSD(track)

MSD(1) = structFunc_struct(track,'Xf',2);
MSD(2) = structFunc_struct(track,'Yf',2);
MSD(3) = structFunc_struct(track,'Zf',2);