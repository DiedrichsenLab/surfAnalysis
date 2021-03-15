function G=surf_makeLabelGifti(data,varargin)
% function G=surf_makeLabelGifti(data,varargin)
% Builds a Label GIFTI structure from scratch for saving 
%
% INPUT: 
%   data: N x Q data matrix: N:number of nodes, Q: number of data series
%
% VARARGIN: 
%   'anatomicalStruct': Anatomical Structure for header
%                       'CortexLeft','CortexRight','Cerebellum' 
%   'labelNames': Cell array of Names for the labels 
%   'columnNames': Cell array of the Names of different data columns 
%   'labelRGBA':  numLabel x 4 matrix of RGB + alpha (opacacity) values
%
% OUTPUT: 
%   G: Gifti object - save with save(G,'filename')
% ---------------------------------------------------------------------------------
anatomicalStruct = 'CortexLeft'; 
columnNames     = {}; 
labelRGBA       = []; 
labelNames      = {}; 
vararginoptions(varargin,{'anatomicalStruct','columnNames','labelNames','labelRGBA'}); 

[N,Q] = size(data);
keys   = unique(data)'; 
nLabel = length(keys); % number of unique labels
% Create naming and coloring if not given:
% 1) Make column_names if empty 
if (isempty(columnNames))
    for i=1:nLabel
        columnNames{i}=sprintf('col_%d',i);
    end 
end 

% 2) Determine color scale if empty
if (isempty(labelRGBA))
    col = hsv(nLabel);
    col = col(randperm(nLabel),:); % shuffle the order so it's more visible
    labelRGBA = zeros(nLabel,4);
    for i=1:nLabel
        labelRGBA(i,:)=[col(i,:) 1];
    end 
end

% 3) Give label names
if (isempty(labelNames))
    for i=1:nLabel
        labelNames{i}=sprintf('label-%d',i);
    end
end


% Create the label.gii structure
this.metadata(1) = struct('name','AnatomicalStructurePrimary','value',anatomicalStruct);  
this.metadata(2) = struct('name','encoding','value','XML_BASE64_GZIP'); 
this.label.name  = labelNames; 
this.label.key   = keys; 
this.label.rgba  = labelRGBA; 
for i=1:Q 
    this.data{i}.data=int32(data(:,i)); 
    this.data{i}.metadata(1) = struct('name','Name','value',columnNames{i});
    this.data{i}.attributes.Dim=N;
    this.data{i}.attributes.DataType = 'NIFTI_TYPE_INT32'; 
    this.data{i}.attributes.Intent   = 'NIFTI_INTENT_LABEL'; 
    this.data{i}.space=[]; 
end; 
G=gifti(this); 