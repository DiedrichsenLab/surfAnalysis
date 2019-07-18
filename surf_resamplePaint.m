function surf_resamplePaint(hemisphere,inFile,outFile,atlasDir,varargin)
% function G=surf_resamplePaint(data,varargin)
% Builds a label.gii structure and resamples it to target
% atlas space.
% / / / / /
% INPUT: 
%   hemisphere :  left(1) or right(2). Currently no support to do both in one call
%   inFile     :  Filename of paint file to be resampled
%   outFile    :  Filename of label.gii file that has been resampled.
%   atlasDir   :  Directory with files for the standard_mesh (Atlas_templates/standard_mesh)
% / / / / /
% VARARGIN: 
%   'resolution' :  Number of vertices in surface which we use to resample
%                  label file. '32k', or '164k' (default = '32k')
% SArbuckle 07/2019

resolution = '32k';
vararginoptions(varargin,{'resolution'});


% error handling
[~,~,ext] = fileparts(inFile);
if ~strcmp(ext,'.paint')
    error('inFile is not a .paint file.')
end
[~,~,ext] = fileparts(outFile);
if~strcmp(ext,'.gii')
    error('outFile is not a .gii file.')
end

% load .metric file 
C = caret_load(inFile);
%nLabel = length(unique(C.data)); % number of unique labels
nLabel = numel(C.paintnames);
% check for invalid text characters in paintnames field
for i = 1:nLabel
    tf = isstrprop(C.paintnames{i},'punct');
    C.paintnames{i}(tf) = '';
end
% set colors
rgb = hsv(nLabel);
rgb = rgb(randperm(nLabel),:); % shuffle the order so it's more visible
% assume empty paintnames ('???') are not labelled and so color set to
% totally transparent
transparent = double(cellfun(@(x) ~isempty(x),C.paintnames))';
rgba = [rgb, transparent];
% convert paint to label.gii
anatomicalStruct = {'CortexLeft','CortexRight'};
G = surf_makeLabelGifti(C.data,...
    'columnNames',C.column_name,...
    'labelNames',C.paintnames,...
    'labelRGBA',rgba,...
    'anatomicalStruct',anatomicalStruct{hemisphere});
save(G,outFile); % convert metric file to gifti and save temporarily (we overwrite this with resampled metric)

% filenames for inflated surface spheres (we align inSphere to outSphere and then resample inFile)
Hem       = {'L','R'}; 
inSphere  = fullfile(atlasDir,['fs_' Hem{hemisphere}],['fsaverage.' Hem{hemisphere} '.sphere.164k_fs_' Hem{hemisphere} '.surf.gii']);
outSphere = fullfile(atlasDir,'resample_fsaverage',['fs_LR-deformed_to-fsaverage.' Hem{hemisphere} '.sphere.' resolution '_fs_LR.surf.gii']);

% do resampling of metric file
system(['wb_command -label-resample ' outFile ' ' inSphere ' ' outSphere ' ADAP_BARY_AREA ' outFile ' -area-surfs ' inSphere ' ' outSphere]);

end