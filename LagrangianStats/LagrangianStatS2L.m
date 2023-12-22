function [S2L] = LagrangianStatS2L(track)

[S2L(1), ~, ~]= structFunc_struct(track,'Vx',2);
[S2L(2), ~, ~]= structFunc_struct(track,'Vy',2);
[S2L(3), ~, ~]= structFunc_struct(track,'Vz',2);
