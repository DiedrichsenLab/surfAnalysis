function surf_icosahedron(atlasDir,varargin)
% function varargout = surf_icosahedron(atlasDir,varargin)
% makes the label file for surface LR
% INPUT:
%           - atlasDir: directory of group FS_LS
% OUTPUT:
%           - label gifti file with tessels
% VARARGIN:
%           - freq: frequency of sampling (2:2:8)
%           - radius
%           - res: resolution (14 or 42 k)
% note:
% 2 - 42 hexagons
% 4 - 162 hexagons
% 6 - 362 hexagons
% 8 - 642 hexagons
% 10 - 1002 hexagons
% 12 - 1442 hexagons
% ---------------------------------------------------
freq    = 6;
radius  = 100;
res     = 164; % using 164 or 32k FS_LR
thres   = 0.25; % how many vertices need to be included (after medial wall removal)
hem = {'L','R'}; % both hemispheres
hemName = {'CortexLeft','CortexRight'};
vararginoptions(varargin,{'freq','radius','res','thres'});
cd(atlasDir);

h=1;
G = gifti(sprintf('fs_LR.%dk.%s.sphere.surf.gii',res,hem{h}));
[x,y,z,~] = make_icosahedron(freq,radius,1,0,1); % make icosahedron
DT = delaunayTriangulation(x',y',z'); % triangulate again (Matlab framework)
[T,Xb] = freeBoundary(DT);
TR = triangulation(T,Xb);
nCentr = size(TR.Points,1);
for i=1:nCentr
    distMat(i,:)=sqrt(sum(bsxfun(@minus,G.vertices,TR.Points(i,:)).^2,2)); % squared Euclidean distance
end
% find minimal distance and assign the label of that centroid
[~,centr]=min(distMat);
% number of vertices per centroid
nVert1 = [[0;unique(centr)'] [0;hist(centr,unique(centr))']]; % add 0
% add a medial wall (overwrite any previous tessel)
M = gifti(sprintf('Yeo_JNeurophysiol11_7Networks.%dk.%s.label.gii',res,hem{h}));
centr(M.cdata==0)=0;
% find now new number of vertices per centroid
nVert2 = [unique(centr)' hist(centr,unique(centr))'];
% remove all other small where less than threshold proportion of
% vertices from the original number
subVert1 = nVert1(ismember(nVert1(:,1),nVert2(:,1)),:); % overlapping rows (centroids still present)
clust = nVert2(:,2)./subVert1(:,2) < thres;
centr(ismember(centr,nVert2(clust,1)))=0;
nLabel = length(unique(centr))-1; % without medial wall
% provide additional info (colour, label, column)
% make the label number continuous -  1:nLabel
% important in case the medial wall wipes out some parcels completely
centr_clean = centr; % clean version - without medial wall
origLab = unique(centr); % original label
% initialise for gifti structure
columnNames     = cell(nLabel+1,1);
labelNames      = cell(nLabel+1,1);
columnNames{1} = 'medial-wall';
labelNames{1}  = 'medial-wall';
labelRGBA       = zeros(nLabel+1,4);
labelRGBA(1,:)  = [0 0 0 1];
col = hsv(nLabel);
col = col(randperm(nLabel),:); % shuffle the order so it's more visible
for i=1:nLabel
    columnNames{i+1}=sprintf('col_%d',i);
    labelRGBA(i+1,:)=[col(i,:) 1];
    labelNames{i+1}=sprintf('label-%d',i);
    centr_clean(centr==origLab(i+1))=i;
end;
for h=1:2
    % create a gifti
    G2 = surf_makeLabelGifti(centr_clean','anatomicalStruct',hemName{h},'columnNames',columnNames,'labelNames',labelNames,'labelRGBA',labelRGBA);
    save(G2,fullfile(atlasDir,sprintf('Icosahedron-%d.%dk.%s.label.gii',nCentr,res,hem{h})));
end

