function surf_groupGiftis(P,varargin)
% function surf_groupGiftis(P,varargin)
% Takes N func.gii files with P columns each (i.e., from each subject) 
% and creates P group func.gii files, with N columns (i.e., for each condition)
% Also saves group summary file (mean across columns), if desired 
% INPUTS: 
%   Infiles: cell array of input files, prompted if empty
% VARARGIN_OPTIONS: 
%   'inputcol',colidx: Only uses the column[colidx] from the input files (DEFAULT: do all)
%   'outfilenames',{'f2','f1'}  :  Name the outputfiles (DAFAULT: uses outfilenamePattern to generate them)
%   'outcolnames',{'n21','n1'}  : Name of columns in output files (DEFAULT: names of input files)
%   'outfilenamePattern',pat    : For automatic naming of group files, it uses sprintf('%s.func.gii',column); 
%   'groupsummary', name : Name of the group summary file that holds each
%   'groupstats', fcn    : Function to the called for group summary file (DEFAULT: nanmean(x,2))
%   'replaceNaNs': NaNs will be replaced with 0 
% joern.diedrichsen@googlemail.com

inputcol=[];
outcolnames={};
outfilenames={};
outfilenamePattern='%s.func.gii';
groupsummary=[];
groupstats = @(x) nanmean(x,2);  

replaceNaNs=0; 
vararginoptions(varargin,{'inputcol','outcolnames','outfilenames','outfilenamePattern','replaceNaNs','groupsummary','groupstats'},{});

if (nargin<1 | isempty(P))
    P=spm_get(inf,{'*.gii'},{'Choose gifti files to summarize'});
end;
if (iscell(P))
    P=char(P);
end;

[a,b,type]=fileparts(P(1,:));
for i=1:size(P,1)
    filename=deblank(P(i,:)); 
    if ~exist(filename,'file') 
        warning(sprintf('File %s does not exist: replaceing with NaNs',filename)); 
        numCols(i)=NaN;
        numRows(i)=NaN;
    else
        M{i}=gifti(filename);
        [numRows(i),numCols(i)]=size(M{i}.cdata);
    end;
end;

% Check if metric files are same format 
if var(numCols(~isnan(numCols)))>0
    warning('Number of columns is not the same');
end;
if var(numRows(~isnan(numRows)))>0
    error('Number of vertices are not the same');
end;

% Determine number of columns 
if (isempty(inputcol))
    inputcol=[1:min(numCols)]; 
end; 

numrows=nanmean(numRows); 

% Get the inputcolumn names 
inputColumnNames = surf_getGiftiColumnNames(M{1}); 

% Get main anatomical structure 
anatomicalStruct = surf_getGiftiAnatomicalStruct(M{1}); 

% Determine output column names 
if (isempty(outcolnames))
    for i=1:length(M)
        [a,b,c]=fileparts(P(i,:));
        outcolnames{i}=b;
    end; 
end; 

% Reorder into new metric files 
for file=1:length(inputcol) 
    OutData=zeros(numrows,length(M))*NaN; 
    for col=1:length(M)
        if (~isempty(M{col}.cdata))           % data is present 
            OutData(:,col)=M{col}.cdata(:,inputcol(file));
        end;
    end;    
    if (isempty(outfilenames) || length(outfilenames)<file) 
        outfilenames{file}=sprintf(outfilenamePattern,inputColumnNames{inputcol(file)}); 
    end;
    O=surf_makeFuncGifti(OutData,'columnNames',outcolnames,'anatomicalStruct',anatomicalStruct);
    save(O,outfilenames{file}); 
    SumData(:,file)=groupstats(OutData); 
end;

% If Summary file name is given, 
if (~isempty(groupsummary))
    O=surf_makeFuncGifti(SumData,'columnNames',inputColumnNames,'anatomicalStruct',anatomicalStruct); 
    save(O,groupsummary); 
end; 

