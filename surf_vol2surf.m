function [G,D]=surf_vol2surf(c1,c2,V,varargin)
% function G=surf_vol2surf(c1,c2,V,varargin)
% Maps functional volume data onto a surface, defined by white and pial
% surface
% INPUTS:
%   c1: Px3 Matrix of Vertices x,y,z coordinates on the white surface
%   c2: Px3 Matrix of Vertices x,y,z coordinates on the pial surface
%   V:  List of volumes (spm_vol) to be mapped
% VARARGIN:
%   'ignore_zeros'       : should zeros be ignored? (default false. set to 1 for F and accuracy)
%   'column_names'       : cell array of column_names for metric file
%   'depths'             : Depths of points along line at which to map 0:white/gray 1: pial
%   'hold'               : Interpolation: 0: nearest neighbour, 1:trilinear
%   'stats'              : Statistics to be evaluated (default: @(x)nanmean(x,2))
%                        : @(x)nanmean(x,2) default and used for activation data 
%                        : @(x)mode(x,2) used when discrete labels are sampled. The most frequent label is assigned 
%   'exclude_thres'      : Threshold enables the exclusion of voxels that touch the surface in two distinct places - 
%                        i.e. Voxels that lie in the middle of a sulcus. If a voxel projects to two separate place 
%                        on the surface, the algorithm excludes it, if the proportion of the bigger cluster
%                        is smaller than the threshold. (i.e. threshold = 0.9 means that the voxel has to
%                        lie at least to 90% on one side of the sulcus) 
%    'faces'             : For threshold exclusion, you need to provide the
%                           faces data from the surface (numFaces x 3 matrix)
%    'anatomicalStruct'  : 'Cerebellum','CortexLeft','CortexRight'
% OUTPUT:
%    M                   : Gifti object- can be saved as a *.func.gii or *.label.gii file 
%    D:         Data structure that contains sparse double [numVox x numNodes] that counts the number of
%               times each voxel is sampled to each node. For example,...
%               ...if voxel 1 is sampled 6x to node 1 (i.e. the same voxel
%               is sampled for each sampling depth), vox2Node(1,1) = 6.
%               Included in this datastructure are the node IDs and the
%               linear voxel indicies. Voxels that are not sampled into
%               surface are dropped from this matrix.
% 
%  Function enables mapping of volume-based data onto the vertices of a
%  surface. For each vertex, the function samples the volume along the line 
%  connecting the white and gray matter surfaces. The points along the line
%  are specified in the variable 'depths'. default is to sample at 5
%  locations between white an gray matter surface. Set 'depths' to 0 to
%  sample only along the white matter surface, and to 0.5 to sample along
%  the mid-gray surface. 
%  The averaging across the sampled points for each vertex is dictated by
%  the variable 'stats'. For functional activation, use 'mean' or
%  'nanmean'. For discrete label data, use 'mode'. 
%  If 'exclude_thres' is set to a value >0, the function will exclude voxels that 
%  touch the surface at multiple locations - i.e. voxels within a sulcus
%  that touch both banks. Set this option, if you strongly want to prevent
%  spill-over of activation across sulci. Not recommended for voxels sizes
%  larger than 3mm, as it leads to exclusion of much data. 
% 
%  For alternative functionality see wb_command volumne-to-surface-mapping 
%  https://www.humanconnectome.org/software/workbench-command/-volume-to-surface-mapping
% 
% 
% joern.diedrichsen@googlemail.com, Feb 2019
% saarbuckle, May 2020: added vox2Node output structure D
ignore_zeros=0;         % In the volume-data, should zero be treated as NaNs? 
exclude_thres=0;        % theshold to exclude voxels that lie within a sulcus 
column_names={};        % Names of the mapped columns (defaults to 
depths=[0 0.2 0.4 0.6 0.8 1];
interp=0;               % Interpolation type 0: Nearest Neighbour, 1: Trilinear
stats=@(x)nanmean(x,2); % What statistics should be used for 
anatomicalStruct='CortexLeft'; 
faces = []; 
vararginoptions(varargin,{'ignore_zeros','column_names','depths','interp','stats','exclude_thres','faces','anatomicalStruct'});

numPoints=length(depths);

nverts=size(c1,2);
if ~isequal([size(c1,2) size(c2,2) ],[3 3])
    error('Coordinates should be Px3 ');
end

if ~isequal(size(c1,2),size(c2,2))
    error('Coordinates C1 and C2 should have same number of elements');
end

% If character array, make into cell
if (ischar(V))
    for i=1:size(V,1);
        A{i}=deblank(V(i,:));
    end;
    V=A;
end;
A=[];
firstgood=[]; 
for i=1:size(V,2);
    try
        A{i}=spm_vol(V{i});
        if (isempty(firstgood))
            firstgood=i;
        end;
    catch
        warning(sprintf('%s could not be opened',V{i}));
        A{i}=[]; 
    end;
end;
V=A;

if (isempty(firstgood))
    error('none of the images could be opened');
end; 


if (length(ignore_zeros)==1)
    ignore_zeros=ones(length(V),1)*ignore_zeros;
end;

for i=1:numPoints % Number of sample points between white and pial surface
    c=(1-depths(i))*c1'+depths(i)*c2';
    indices(:,i)=surfing_coords2linvoxelidxs(c,V{firstgood});
end;

% If necessary, now ensure that voxels are mapped to on continuous location
% only on the flat map - exclude voxels that are assigned to two sides of
% the sulcus
if (exclude_thres>0)
    exclude=zeros(prod(V{1}.dim),1);
    if (isempty(faces))
        error('provide topology data (faces), so that projections should be avoided');
    end;
    S.Tiles.data=double(faces); % GA: leaving without conversion to double returns error with function sparse below
    
    % Precaluclate the Edges for fast cluster finding
    fprintf('Calculating Edges\n');
    S.num_nodes=max(max(S.Tiles.data));
    S.Edges.data=[];
    for i=1:3
        i1=S.Tiles.data(:,i);i2=S.Tiles.data(:,mod(i,3)+1);
        S.Edges.data=[S.Edges.data;[i1 i2]];
    end;
    S.Edges.data=unique(S.Edges.data,'rows');
    
    % Generate connectivity matrix
    G=sparse(S.Edges.data(:,1),S.Edges.data(:,2),ones(size(S.Edges.data,1),1),S.num_nodes,S.num_nodes);
    
    % Cluster the projections to the surface and exclude all voxels
    % in which the weight of the biggest cluster is not > thres
    fprintf('Checking projections\n');
    I=unique(indices(~isnan(indices)))';
    for i=I
        
        % Calculate the weight of voxel on node
        weight=sum(indices==i,2);
        indx=find(weight>0);
        
        % Check whether the nodes cluster
        H=G(indx,indx);
        H=H+H'+speye(size(H,1));
        [p,q,r,s]=dmperm(H);
        CL=zeros(size(indx));
        for c=1:length(r)-1
            CL(p(r(c):r(c+1)-1),1)=c;
        end;
        
        
        if (max(CL)>1)
            weight_cl=zeros(max(CL),1);
            for cl=1:max(CL)
                weight_cl(cl,1)=sum(weight(indx(CL==cl)));
            end;
            [m,cl]=max(weight_cl);
            if (m/sum(weight_cl)>exclude_thres)  % Can be assigned to one of the clusters: exclude from others
                A=indices(indx(CL~=cl),:);
                A(A==i)=NaN;
                indices(indx(CL~=cl),:)=A;
                exclude(i)=1;
                % fprintf('assigned: %2.3f\n',m/sum(weight_cl));
            else                                 % Cannot be assigned - kill completely
                A=indices(indx,:);
                A(A==i)=NaN;
                indices(indx,:)=A;
                exclude(i)=2;
                % fprintf('excluded: %2.3f %d \n',m/sum(weight_cl),max(CL));
            end;
        end;
    end;
    
    % For debugging: save the volume showing excluded voxels in current
    % directory
    Vexcl=V{1};
    Vexcl.fname='excl.nii';
    spm_write_vol(Vexcl,reshape(exclude,Vexcl.dim));
end;

i=find(~isnan(indices));
%data=zeros(size(indices))*NaN; %GA: Workbench doesn't know to ignore NaNs
%(as SPM does), but it ignores zeros instead (see changes below)
data=zeros(size(indices));
for v=1:length(V);
    if (~isempty(V{v}))
        X=spm_read_vols(V{v});
        %         if (ignore_zeros(v)) %GA
        %             X(X==0)=NaN;
        %         end;
        X(isnan(X))=0; %GA
        data(i)=X(indices(i)); %GA
        M.data(:,v)=stats(data);
    else 
        M.data(:,v)=nan(size(c1,1),1); 
    end; 
end;
fprintf('\n');
% Calc vox2node sparse matrix: (SA 05/2020)
% We will find which node each voxel maps onto, re-counting a voxel if it
% is mapped multiple times to same node (b/c sampling can occur at multiple
% depths)
maxVoxID  = numel(X);   % # unique voxels
maxNodeID = size(c1,1); % # unique nodes
matNodes  = repmat(1:maxNodeID,numPoints,1)';   % vector of node ids sized to match indices matrix
vox2Node  = sparse(double(indices(i)), matNodes(i), ones(size(i)), maxVoxID, maxNodeID); 
vox2Node  = vox2Node * diag(sparse(1./sum(vox2Node,1))); % normalize each column (node) to 1:
% Some nodes are empty, but don't drop those from the metric file.
% However, we drop empty voxels from the sparse connection matrix:
keepVox = find(sum(vox2Node,2)~=0);
D.vox2Node = vox2Node(keepVox,:);
D.linVoxID = keepVox';
D.nodeID   = 1:maxNodeID;
% vox2Node'*X(:); % to check. BUT, not that vox2Node cannot deal with
% voxels with vals==nan

% Determine the column names based on the filenames of the volumes
if (isempty(column_names))
    for i=1:length(V)
        if (~isempty(V{i}))
            [~,column_names{i},~]=spm_fileparts(V{i}.fname);
        else 
            column_names{i}='void';
        end; 
    end;
end;

% Build a functional GIFTI structure
this.metadata(1) = struct('name','AnatomicalStructurePrimary','value',anatomicalStruct);  
this.metadata(2) = struct('name','encoding','value','XML_BASE64_GZIP'); 
this.label.name  = {'???'}; 
this.label.key   = 0; 
this.label.rgba  = [1 1 1 0]; 
for i=1:length(V) 
    this.data{i}.data=single(M.data(:,i)); 
    this.data{i}.metadata(1) = struct('name','Name','value',column_names{i});
    this.data{i}.attributes.Dim=size(M.data,1);
    this.data{i}.attributes.DataType = 'NIFTI_TYPE_FLOAT32'; 
    this.data{i}.attributes.Intent   = 'NIFTI_INTENT_NONE'; 
    this.data{i}.space=[]; 
end; 
G=gifti(this); 