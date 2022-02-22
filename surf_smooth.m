function Out=surf_smooth(P,varargin)
% Smooth metric files 
% Out=surf_smooth(P,varargin);
% Uses wb_command to smooth the files 
% INPUT: 
%   P: cell array of metric files to smooth 
% VARARGIN: 
%   'surf',file:    Surface file to smooth on (fyll filepath)
%   'kernel':       The sigma for the gaussian kernel function, in mm
% OUTPUT: 
%   Cell array of smooth metric files 
% -------------------------------------------------------------------
% Note: default algorithm GEO_GAUSS_AREA, other methods not implemented
% More info: https://www.humanconnectome.org/software/workbench-command/-metric-smoothing
% EBerlot, May 2019
%
prefix  = 's';
kernel  = 2.0;
surf    = [];
vararginoptions(varargin,{'prefix','kernel','surf','outfile'});
if (nargin<1 || isempty(P))
    error('Provide the metric .func.gii file.\n');
elseif (iscell(P))
    P=char(P);
end;
if (isempty(surf))
    error('Provide the .surf.gii file.\n');
end;
N   = size(P,1);
Out = cell(1,N);
for i=1:N 
    [dir,name,postfix]=fileparts(deblank(P(i,:)));
    Out{i}=fullfile(dir,[prefix name postfix]);
    comm=sprintf('wb_command -metric-smoothing %s %s %f %s -fix-zeros',...
        surf,P,kernel,Out{i});
    %fprintf('%s\n',comm) 
    [err,out]=system(comm);
    fprintf('\nSmoothed file %s\n%d:%s\n',[name postfix],err,out);
    %fprintf('Smoothed file %s\n',[name postfix]);
end;

