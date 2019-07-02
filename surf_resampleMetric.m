function surf_resampleMetric(hemisphere,inFile,outFile,atlasDir,varargin)
% function G=surf_resampleMetric(data,varargin)
% Builds a functional funtional GIFTI structure and resamples it to target
% atlas space.
% / / / / /
% INPUT: 
%   hemisphere :  left(1) or right(2). Currently no support to do both in one call
%   inFile     :  Filename of metric file to be resampled
%   outFile    :  Filename of metric file that has been resampled.
%   atlasDir   :  Directory with files for the standard_mesh (Atlas_templates/standard_mesh)
% / / / / /
% VARARGIN: 
%   'surf'     :  Number of vertices in surface which we use to resample
%                  metric file. '32k', or '164k' (default = '32k')
% SArbuckle 05/2019

resolution = '32k';
vararginoptions(varargin,{'resolution'});


% error handling
[~,~,ext] = fileparts(inFile);
if ~strcmp(ext,'.metric')
    error('inFile is not a .metric file.')
end
[~,~,ext] = fileparts(outFile);
if~strcmp(ext,'.gii')
    error('outFile is not a .gii file.')
end

% load .metric file and convert to gifti
anatomicalStruct = {'CortexLeft','CortexRight'};
C = caret_load(inFile);
G = surf_makeFuncGifti(C.data,'columnNames',C.column_name,'anatomicalStruct',anatomicalStruct{hemisphere});
save(G,outFile); % convert metric file to gifti and save temporarily (we overwrite this with resampled metric)

% filenames for inflated surface spheres (we align inSphere to outSphere and then resample inFile)
Hem       = {'L','R'}; 
inSphere  = fullfile(atlasDir,['fs_' Hem{hemisphere}],['fsaverage.' Hem{hemisphere} '.sphere.164k_fs_' Hem{hemisphere} '.surf.gii']);
outSphere = fullfile(atlasDir,'resample_fsaverage',['fs_LR-deformed_to-fsaverage.' Hem{hemisphere} '.sphere.' resolution '_fs_LR.surf.gii']);

% do resampling of metric file
system(['wb_command -metric-resample ' outFile ' ' inSphere ' ' outSphere ' ADAP_BARY_AREA ' outFile ' -area-surfs ' inSphere ' ' outSphere]);

end