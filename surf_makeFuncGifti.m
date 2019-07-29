function G=surf_makeFuncGifti(data,varargin)
% function G=surf_makeFuncGifti(data,varargin)
% Builds a functional GIFTI structure from scratch for saving 
% INPUT: 
%   data: N x Q data matrix: N:number of nodes, Q: number of data series
% VARARGIN: 
%   'anatomicalStruct': Anatomical Structure for header
%                       'CortexLeft','CortexRight','Cerebellum' 
%   'columnNames':    : Cell array of the names of different columns 
% OUTPUT: 
%   G: Gifti object - save with save(G,'filename')

anatomicalStruct = 'CortexLeft'; 
columnNames     = {}; 

vararginoptions(varargin,{'anatomicalStruct','columnNames'}); 

[N,Q] = size(data);
% Make column_names if empty 
if (isempty(columnNames))
    for i=1:Q 
        columnNames{i}=sprintf('col_%d',i);
    end; 
end; 
        
this.metadata(1) = struct('name','AnatomicalStructurePrimary','value',anatomicalStruct);  
this.metadata(2) = struct('name','encoding','value','XML_BASE64_GZIP'); 
this.label.name  = {'???'}; 
this.label.key   = 0; 
this.label.rgba  = [1 1 1 0]; 
for i=1:Q 
    this.data{i}.data=single(data(:,i)); 
    this.data{i}.metadata(1) = struct('name','Name','value',columnNames{i});
    this.data{i}.attributes.Dim=N;
    this.data{i}.attributes.DataType = 'NIFTI_TYPE_FLOAT32'; 
    this.data{i}.attributes.Intent   = 'NIFTI_INTENT_NONE'; 
    this.data{i}.space=[]; 
end; 
G=gifti(this); 