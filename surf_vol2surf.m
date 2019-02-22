function M=surf_vol2surf(c1,c2,V,varargin)
% function M=caret_vol2surf_own(c1,c2,V,varargin)
% Maps functional volume data onto a surface, defined by white and pial
% surface
% INPUTS:
%   c1: Px3 Matrix of Node x,y,z coordinates on the white surface
%   c2: Px3 Matrix of Node x,y,z coordinates on the pial surface
%   V: List of volumes (spm_vol) to be mapped
% VARARGIN:
%   'ignore_zeros',0     : default none: should zeros be ignored (set to 1 for F and accuracy)
%   'column_names'       : cell array of column_names for metric file
%   'depths'             : Depth on the line for which to map 0:white/gray 1: pial
%   'hold'               : Interpolation: 0: nearest neighbour, 1:trilinear
%   'stats'              : Statistics to be evaluated (default: @(x)nanmean(x,2))
%                        : @(x)nanmean(x,2) default and used for activation data 
%                        : @(x)mode(x,2) used when discrete labels are sampled. The most frequent label is assigned 
%   'surf',surf          : Precaluclated standard surface for clustering
%   'exclude_thres',thres: Threshold enables the exclusion of voxels that touch the surface in two distinct places - 
%                        i.e. Voxels that lie in the middle of a sulcus. The threshold is the proportion of the bigger cluster
%                        necessary. (i.e. 0.9 means that the voxel has to
%                        lie at least to 90% on one side of the sulcus to be included in the mapping
%                        0 maps all voxels. 
% OUTPUT:
%    M                   : Gifti object- can be saved as a *.func.gii or *.label.gii file 
% joern.diedrichsen@googlemail.com, Feb 2019

ignore_zeros=0;         % In the volume-data, should zero be treated as NaNs? 
exclude_thres=0;        % theshold to exclude voxels that lie within a sulcus 
column_names={};        % Names of the mapped columns (defaults to 
depths=[0 0.2 0.4 0.6 0.8 1];
interp=0;               % Interpolation type 0: Nearest Neighbour, 1: Trilinear
stats=@(x)nanmean(x,2); % What statistics should be used for 

vararginoptions(varargin,{'ignore_zeros','column_names','depths','interp','stats','exclude_thres','topo'});

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
    if (isempty(topo))
        error('provide topo file, if dual projections should be avoided');
    end;
    if ischar(topo)
        topo=caret_load(topo);
        S.Tiles.data=topo.data;
    else
        S.Tiles.data=topo;
    end;
    
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
data=zeros(size(indices))*NaN;

for v=1:length(V);
    if (~isempty(V{v}))
        X=spm_read_vols(V{v});
        if (ignore_zeros(v))
            X(X==0)=NaN;
        end;
        data(i)=X(indices(i));
        M.data(:,v)=stats(data);
    else 
        M.data(:,v)=nan(size(c1,1),1); 
    end; 
end;

% Determine the column names based on the filenames of the volumes
if (isempty(column_names))
    for i=1:length(V)
        if (~isempty(V{i}))
            [~,column_names{i},~]=spm_fileparts(V{i}.fname);
        else 
            column_name{i}='void';
        end; 
    end;
end;
M=caret_struct('metric','data',M.data,'column_name',column_names);
