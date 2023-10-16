%% ----------------------------- DESCRIPTION ----------------------------- &&
% Convert EBSD data in ctf format into a valid DAMASK simulation domain

%% ----------------------------- STARTUP AND IMPORT -----------------------%%
addpath('..\..\..\01_SOFTWARE_ENGINE_ROOM\00_MATLAB\00_MTEX\mtex-5.9.0\mtex-5.9.0')
startup_mtex

disp('Reading ctf file ...')
CS = {... 
  'notIndexed',...
  'notIndexed',...
  'notIndexed',...
  'notIndexed',...
  'notIndexed',...
  crystalSymmetry('6/mmm', [3 3 4.7], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0 0 0.55]),...
  'notIndexed'};

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

fname = '..\01_geom\CTF_File_68.ctf';
ebsd = EBSD.load(fname,CS,'interface','ctf');
ebsd = ebsd('Ti-Hex');

%% ----------------------------- CROPPING ---------------------------------%%
disp('Cropping EBSD image ...')

x1 = 45;
x2 = 72;
y1 = 0;
y2 = 20;

ebsd = ebsd((x1<ebsd.x) & (ebsd.x<x2) & (y1<ebsd.y) & (ebsd.y<y2));
plot(ebsd,ebsd.orientations)
%% ----------------------------- GRAIN RECONSTRUCTION ---------------------%%
disp('Reconstructing grains ...')
alpha          = 5;
size_threshold = 5;
[grains,ebsd.grainId,ebsd.mis2mean]              = calcGrains(ebsd,'angle',alpha*degree);
% plot(grains,grains.meanOrientation)

disp('Remove small grains ...')
ebsd(grains(grains.grainSize<=size_threshold))   = [];
[grains,ebsd.grainId,ebsd.mis2mean]              = calcGrains(ebsd,'angle',alpha*degree);
% plot(grains,grains.meanOrientation)
%% ----------------------------- FILL EBSD --------------------------------%%
disp('Fill EBSD ...')

ebsd=fill(ebsd);
grain_Id = ebsd.grainId;

% Crop edges as long as there are pixels without grain ID
idn=any(isnan(grain_Id(:)));
while idn
    grain_Id([1,end],:)=[];
    grain_Id(:,[1,end])=[];
    idn=any(isnan(grain_Id(:)));
end


%% ----------------------------- EXPORT TO DAMASK SYNTAX ------------------%%
disp('Export data to Damask syntax ...')
filename = '../01_geom/grid_grain_IDs.txt';
writematrix(transpose(grain_Id),filename)

eulers = rad2deg([grains.meanOrientation.phi1 grains.meanOrientation.Phi grains.meanOrientation.phi2]);

filename = '../01_geom/grain_euler_angles.txt';
writematrix(round(eulers,2),filename,Delimiter=',')

filename = '../01_geom/dimensions.txt';
writematrix(round([ebsd.xmin ebsd.xmax ebsd.ymin ebsd.ymax ebsd.dx ebsd.dy],4),filename,Delimiter=',')
