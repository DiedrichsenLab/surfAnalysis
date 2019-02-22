function S=surf_resliceFS2WB(subj_name,subject_dir,atlas_dir,out_dir,varargin);
% functionsurf_resliceFS2WB(subj,subject_dir,atlas_dir,out_dir,varargin);
% Resampels a registered subject surface from freesurfer average to the new
% symmetric fs_LR_164 surface, standard in workbench.  
% This allows things to happen exactly in atlas space - each vertex number
% corresponds exactly to a anatomical location 
% INPUT: 
%   subj: subject name
%   subjects_dir: freesurfer's SUBJECT_DIR
%   atlas_dir: Directory with files for the standard_mesh (Atlas_templates/standard_mesh)
%   out_dir: The new resampled files will be placed in (out_dir/subj)
% VARARGIN: 
%   'hemisphere',[1 2]  : left / right or both hemispheres 
%   'align_surf':       : Shift the surface to correct for freesurfer convention? 
%   'surf_files',{'',''}: Surface files to be resampled 
%                       {'.white','.pial','.inflated','.sphere.reg','.sphere'};
%   'curv_files',{'',''}: Curvature files to be resampled 
%                       {'.curv','.sulc','.area'}
% ---------------------------
current_dir=pwd;

direct=[getenv('FREESURFER_HOME') filesep 'average' filesep 'surf' ];
if (isempty(subject_dir))
    subject_dir=getenv('SUBJECTS_DIR'); 
end; 
structname={'left','right'};
hem={'lh','rh'};
Hem={'L','R'}; 
hemisphere=[1:2]; % Do both hemispheres 
surf_files={'.white','.pial','.inflated'};
curv_files={'.curv','.sulc','.area'}; 
smoothing=1; 
align_surf=[1 1 1]; 

vararginoptions(varargin,{'smoothing','surf_files','curv_files','hemisphere','align_surf'});

% ----------------------------------------------------
% read freesurfer version
ver_file = fullfile(getenv('FREESURFER_HOME'),'build-stamp.txt');
if exist(ver_file)
    fid     = fopen(ver_file);
    verStr  = fgetl(fid);
    fclose(fid);
    verStr  = strrep(verStr,'-v',' ');
    verStr  = textscan(verStr,'%s%f');
    fsl_ver = verStr{2};
end

% ------------------------------------------------------ 
% create new output directory 
new_dir = [out_dir filesep subj_name];
if (~exist(new_dir))
    mkdir(new_dir);
end;

% -----------------------------------------------------
% Transform surface files to standard mesh 
num_surf=length(surf_files); 
cd(fullfile(subject_dir,subj_name,'surf')); 

% -----------------------------------------------------
% Figure out the shifting of coordinate systems:
% Freesurfer uses vertex coordinates in respect to
% the center of the 256x256x256 image.
% Independent of the real zero point in the original image
% So to find a transform of the
% Mvox2surf: Transform of voxels in 256x256 image to surface vertices
% Mvox2space: Transform of voxel to subject space
anafile=fullfile(subject_dir,subj_name,'mri','brain.mgz');
[status,result] = system(['mri_info ' anafile ' --vox2ras-tkr']);
A=sscanf(result,'%f');
Mvox2surf=reshape(A,4,4)';  
[status,result] = system(['mri_info ' anafile ' --vox2ras']);
A=sscanf(result,'%f');
Mvox2space=reshape(A,4,4)';
Msurf2space=Mvox2space*inv(Mvox2surf);
    

% -----------------------------------------------------
% Transform the surfaces from the two hemispheres 
for h=hemisphere
    % Convert reg-sphere 
    reg_sphere = [hem{h} '.sphere.reg.surf.gii']; 
    system(['mris_convert ' hem{h} '.sphere.reg ' reg_sphere]);
    
    % -----------------------------------------------------
    % Do all the surface files 
    for i=1:length(surf_files)    
        % Set up file names 
        file_name = [hem{h} surf_files{i} '.surf.gii']; 
        out_name = fullfile(new_dir,[subj_name '.' Hem{h} surf_files{i} '.164k.surf.gii']); 
        atlas_name = fullfile(atlas_dir,'resample_fsaverage',['fs_LR-deformed_to-fsaverage.' Hem{h} '.sphere.164k_fs_LR.surf.gii']);
        
        % Convert surface to Gifti 
        system(['mris_convert ' hem{h} surf_files{i} ' ' file_name]);
        
        system(['wb_command -surface-resample ' file_name ' ' reg_sphere ' ' atlas_name ' BARYCENTRIC ' out_name]);

        A=gifti(out_name); 
        if (align_surf(i)) 
            [A.vertices(:,1),A.vertices(:,2),A.vertices(:,3)]=spmj_affine_transform(A.vertices(:,1),A.vertices(:,2),A.vertices(:,3),Msurf2space);
        end; 
        save(A,out_name);
    end; 
    
    % -----------------------------------------------------
    % Do all the curvature files 
    for i=1:length(curv_files)    
        % Set up file names 
        file_name = [hem{h} curv_files{i} '.shape.gii']; 
        out_name = fullfile(new_dir,[subj_name '.' Hem{h} curv_files{i} '.164k.shape.gii']); 
        atlas_name = fullfile(atlas_dir,'resample_fsaverage',['fs_LR-deformed_to-fsaverage.' Hem{h} '.sphere.164k_fs_LR.surf.gii']);
        
        % Convert surface to Gifti 
        system(['mris_convert -c ' [hem{h} curv_files{i}] ' ' [hem{h} surf_files{1}] ' ' file_name]);
        system(['wb_command -metric-resample ' file_name ' ' reg_sphere ' ' atlas_name ' BARYCENTRIC ' out_name]);
    end; 
    
end; 


