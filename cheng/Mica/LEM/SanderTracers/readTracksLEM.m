function part = readTracksLEM(filename,fileid)
% 
% part = readTracksLEM(filename)
%
% reads .bin files tracked by Sander in the LEM
% 
% MB 01/02/2018


fid = fopen(filename,'r');
skip_trackid = 16;
skip_time = 18 ;
skip_matchid = 18;
skip_xyz = 16;

fseek(fid,0,-1);
trackid = fread(fid,'uint32',skip_trackid);

fseek(fid,4,-1);
time = fread(fid,'uint16',skip_time);

fseek(fid,6,-1);
matchid = fread(fid,'uint16',skip_matchid);

fseek(fid,8,-1);
x = fread(fid,'float',skip_xyz);

fseek(fid,12,-1);
y = fread(fid,'float',skip_xyz);

fseek(fid,16,-1);
z = fread(fid,'float',skip_xyz);

fclose(fid);

part.file = filename ;
part.matchid = matchid ;
part.t = time ;
part.Ntrack = trackid + i*fileid ;
part.x = x ;
part.y = y ;
part.z = z ;

% %%
% 
% pa
% 
% 
%    BinaryReadList[
%     fn, {"UnsignedInteger32", "UnsignedInteger16", 
%      "UnsignedInteger16", "Real32", "Real32", "Real32"}]; (* trackid, 
%   time, match id, x,y,z,*)
%   
%   tracks = SortBy[GatherBy[tracks, First], {First@*First, Length}];
%   tracks = Apply[{#1, #2, {#4, #5, #6}} &, tracks, {2}]; (* {trackid, 
%   time, {x,y,z}}*)
%   tracks
 