function [Y, P, coord] = surf_cross_section(surface, data, varargin)
%% function [Y, P, coord] = surf_cross_section(surface, data, varargin);
% Samples data as cross-section from either a .border file or a virtual
% strip (done in WorkBench).
%
% INPUT:
%   surface :   Name of surface file where the borders were drawn on. Flat
%               coordinates are (Nx3), with the last column being all
%               zeros.
%   data    :   Name of gifti file to be sampled from. Can be metric,
%               coord, RBG, sulc (nii), or paint like type of file.
%
% VARARGIN:
%   border  :   Name of border file (may contain more than one border). If
%               not specified, it will use the virtual strip option (see
%               below in VARARGIN) with default coordinates instead.
%   b_name  :   If provided (string), only the border specified by b_name
%               is used for cross-section. By default (not specified),
%               every border in border file is used.
%   from,to :   If provided, this function will ignore (nor need) a
%               specified border. It will instead use a virtual strip for
%               the cross-section on *flat.surf.gii defined from flat
%               coordinates "from" and "to" (1x2 array each), with the
%               width of the strip being defined by "width" (in mm) above
%               and below the virtual line where we sample from.
%   width   :   Either width (in mm) of the sampling along the provided
%               border (default is 1 mm), or width of the strip being
%               defined by "from" and "to" coordinates.
%   n_point :   Number of points on for the sampling on the virtual strip
%               (default is 101).
%   outname :   If provided output file name, this function will generate a
%               metric file of the border vertices. By default, the file is
%               generated with the same name as the border file, and
%               automatically deleted.
%   con_map :   Specifies which of the contrast maps (one per subject) to
%               sample from (default is all of the subj maps in .gii file).
%
% OUTPUT:
%   Y       :   sampled data (NxQ) along the borders vertices
%   P       :   logicial value of which vertices are included
%   coord   :   mean flat coordinates of the sampled strip
%
% To get the coordinates for the sampled points, sample also from the
% respective coord file.
% 
% EXAMPLE: 
% surf_cross_section('flat.surf.gii','func.gii',...
%             'from',[-23 80],'to',[124 85], ...
%             'width',10);
% _________________________________________________________________
% Adapted from caret_crosssection and including options from
% caret_crosssection_flat (Joern Diedrichsen 2013)
% Last change: GA - 2020.11.13

% Defaults
border  = [];
outname = [];
con_map = [];
b_name  = [];
from    = [];
to      = [];
width   = 1;
n_point = 101;
vararginoptions(varargin, {'border', 'outname', 'con_map', 'b_name', 'from', 'to', 'width', 'n_point'});

if isempty(surface) || isempty(data)
    error('This function needs a flat surface (.gii) and a metric file (.gii) to sample from.');
end

% Convert niftis into giftis (for sulcal depth info)
[dir,name,ext]=fileparts(data);
if strcmp(ext, '.nii')
    cifti_in  = [dir, '/', name, ext];
    gifti_out = [dir, '/', name '.gii'];
    system(['wb_command -cifti-convert -to-gifti-ext ' cifti_in ' ' gifti_out]);
    data = gifti_out;
end

% Decide whether to use border or virtual strip option
if ~isempty(border) % Border file
    % Decide whether to save border file as func.gii or not (default)
    % -border-to-vertices needs outname, so give same as border file
    if isempty(outname)
        [p,n,~] = fileparts(border);
        outname = fullfile(p, sprintf('%s.func.gii', n));
        saveout = 0;
    else
        saveout = 1;
    end
    
    % Convert border to gifti (sample vertices): outputs a metric with 1s on
    % vertices that follow a border, and 0s elsewhere. By default, a separate
    % metric column is created for each border.
    if ~isempty(b_name)
        system(['wb_command -border-to-vertices ' surface ' ' border ' ' outname ' -border ' b_name]);
    else
        system(['wb_command -border-to-vertices ' surface ' ' border ' ' outname]);
    end
    
    % Load border metric file
    B = gifti(outname);
    
    % Load data metric file
    D = gifti(data);
    
    % Load flat surface coordinate data
    flat_coord = gifti(surface);
    T = flat_coord.vertices;
    
    %     % Perform cross-section
    %     if ~isempty(con_map)
    %         Y = D.cdata(B.cdata==1, con_map);
    %     else
    %         Y = D.cdata(B.cdata==1, :);
    %     end
    %     P = B.cdata;
    %     coord = T(B.cdata==1, :);
    
    % Find points on the border
    idx = find(B.cdata==1);
    n_point = numel(idx);
    
    % Preallocate
    if ~isempty(con_map)
        Y = zeros(n_point-1, numel(con_map));
    else
        Y = zeros(n_point-1, size(D.cdata,2));
    end
    P = zeros(size(T,1),1);
    coord = zeros(n_point-1, size(T,2));
    
    % For every 2 points compute the distance and cross section
    for i=1:n_point-1
        from     = T(idx(i), 1:2);
        to       = T(idx(i+1), 1:2);
        vec      = to-from;  % This is the strip vector
        points   = bsxfun(@minus, T(:,1:2), from);
        project  = (points*vec') ./ (vec*vec');
        residual = points - bsxfun(@times, project, vec);
        distance = sqrt(sum(residual.^2,2));
        
        % Find the vertices on the border
        x           = linspace(0,1,n_point);
        indx        = find(distance<width & project>=x(i) & project<=x(i+1) & sum(T(:,1:2).^2,2)>0);
        % Perform cross-section (contrast map optional)
        if ~isempty(con_map)
            Y(i,:)  = nanmean(D.cdata(indx, con_map));
        else
            Y(i,:)  = nanmean(D.cdata(indx, :));
        end
        P(indx,1)   = 1;
        coord(i,:)  = nanmean(T(indx,:));
    end
    notnan = ~isnan(Y(:,1));
    Y =         Y(notnan, :);
    coord = coord(notnan, :);
    
    % Delete border metric file
    if saveout==0
        delete(outname);
    end
    
else % Virtual strip
    if ~isempty(from)
        if isempty(to)
            error('If you specified the "from" coordinates, you also need to specify the "to" coordinates.');
        end
    else % (Defaults)
        warning('No border file or coordinates provided! Going with default option of virtual strip (10 mm) from coord [-23 80] to coord [124 85].');
        from  = [-23 80];
        to    = [124 85];
        width = 10;
    end
    
    % Load data metric file
    D = gifti(data);
    
    % Load flat surface coordinate data
    flat_coord = gifti(surface);
    T = flat_coord.vertices;
    
    % Find points on the virtual strip
    vec      = to-from;  % This is the strip vector
    points   = bsxfun(@minus, T(:,1:2), from);
    project  = (points*vec') ./ (vec*vec');
    residual = points - bsxfun(@times, project, vec);
    distance = sqrt(sum(residual.^2,2));
    
    % Find the vertices on the virtual strip, while excluding the
    % vertices with (0,0), as these are not mapped onto the flatmap
    x = linspace(0,1,n_point);
    if ~isempty(con_map)
        Y = zeros(n_point-1, numel(con_map));
    else
        Y = zeros(n_point-1, size(D.cdata,2));
    end
    P = zeros(size(T,1),1);
    coord = zeros(n_point-1, size(T,2));
    for i=1:n_point-1
        indx        = find(distance<width & project>=x(i) & project<=x(i+1) & sum(T(:,1:2).^2,2)>0);
        % Perform cross-section (contrast map optional)
        if ~isempty(con_map)
            Y(i,:)  = nanmean(D.cdata(indx, con_map));
        else
            Y(i,:)  = nanmean(D.cdata(indx, :));
        end
        P(indx,1)   = 1;
        coord(i,:)  = mean(T(indx,:));
    end
end


