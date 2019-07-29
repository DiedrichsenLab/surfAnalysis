function G=surf_makeIndexVol(c1,c2,V,filename,varargin)
% function G=surf_makeIndexVol(c1,c2,V,varargin)
%   Maps functional volume data onto a surface, defined by white and pial
%   surface
% INPUTS:
%   c1: Px3 Matrix of Node x,y,z coordinates on the white surface
%   c2: Px3 Matrix of Node x,y,z coordinates on the pial surface
%   V:  Volumne information structure (with filename)
% VARARGIN:
%   'depths'             : Depth on the line for which to map 0:white/gray 1: pial
%   'exclude_thres',thres: Threshold enables the exclusion of voxels that touch the surface in two distinct places -
%                        i.e. Voxels that lie in the middle of a sulcus. The threshold is the proportion of the bigger cluster
%                        necessary. (i.e. 0.9 means that the voxel has to
%                        lie at least to 90% on one side of the sulcus to be included in the mapping
%                        0 maps all voxels.
%    'anatomicalStruct'  : 'Cerebellum','CortexLeft','CortexRight'
% OUTPUT:
%    M                   : Gifti object- can be saved as a *.func.gii or *.label.gii file
% joern.diedrichsen@googlemail.com, Feb 2019

exclude_thres=0;        % theshold to exclude voxels that lie within a sulcus
depths=[0:1/20:1];
vararginoptions(varargin,{'ignore_zeros','column_names','depths','interp','stats','exclude_thres','topo','anatomicalStruct'});

numPoints=length(depths);

nverts=size(c1,1);
if ~isequal([size(c1,2) size(c2,2) ],[3 3])
    error('Coordinates should be Px3 ');
end

if ~isequal(size(c1,1),size(c2,1))
    error('Coordinates C1 and C2 should have same number of elements');
end

% Make a volume of zeros
V.dat = zeros(V.dim);

% Translate all the points into linear voxel indices
for i=1:numPoints % Number of sample points between white and pial surface
    c=(1-depths(i))*c1'+depths(i)*c2';
    indices(:,i)=surfing_coords2linvoxelidxs(c,V);
end;

% If necessary, now ensure that voxels are mapped to on continuous location
% only on the flat map - exclude voxels that are assigned to two sides of
% the sulcus: This needs to be checked still
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
    
    
end;
% Now do the searching
uniqueInd=unique(indices(~isnan(indices))); 
for i=uniqueInd'
    ii=find(indices==i); % Find vertices that have this voxel
    [v,f]=ind2sub(size(indices),ii);
    V.dat(i)=mode(v); % Assign the most frequently occurring index
end;
% For debugging: save the volume showing excluded voxels in current
% directory

V.pinfo=[1;0;0];
V.fname=filename; 
spm_write_vol(V,V.dat);