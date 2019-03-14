function G=surf_makeFuncGifti(data,varargin)
% function G=surf_makeFuncGifti(data,varargin)
% Builds a functional funtional GIFTI structure from scratch for saving 
% INPUT: 
%   data: N x Q data matrix: N:number of nodes, Q: number of data series
% VARARGIN: 
%   'anatomicalStruct': Anatomical Structure for header
%                       'CortexLeft','CortexRight','Cerebellum' 
% 
% OUTPUT: 
%   G: Gifti object - save with save(G,'filename')

anatomicalStruct = 'CortexLeft'; 
column_names     = {}; 

vararginoptions(varargin,{'anatomicalStruct','column_names'}); 

[N,Q] = size(data);
% Make column_names if empty 
if (isempty(column_names))
    for i=1:Q 
        column_names{i}=sprintf('col_%d',i);
    end; 
end; 
        
this.metadata(1) = struct('name','AnatomicalStructurePrimary','value',anatomicalStruct);  
this.metadata(2) = struct('name','encoding','value','XML_BASE64_GZIP'); 
this.label.name  = {'???'}; 
this.label.key   = 0; 
this.label.rgba  = [1 1 1 0]; 
for i=1:Q 
    this.data{i}.data=single(data(:,i)); 
    this.data{i}.metadata(1) = struct('name','Name','value',column_names{i});
    this.data{i}.attributes.Dim=N;
    this.data{i}.attributes.DataType = 'NIFTI_TYPE_FLOAT32'; 
    this.data{i}.attributes.Intent   = 'NIFTI_INTENT_NONE'; 
    this.data{i}.space=[]; 
end; 
G=gifti(this); 