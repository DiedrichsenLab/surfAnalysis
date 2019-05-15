function Y = surf_cross_section(border, surface, data, varargin)
% Samples data as cross-section from a .border file in WorkBench
%
% INPUT:
%   border  :   Name of border file (may contain more than one border).
%   surface :   Name of surface file where the borders were drawn on.
%   data    :   Name of gifti file to be sampled from.
%
% VARARGIN:
%   outname :   If provided output file name, this function will generate a
%               metric file of the border vertices. By default, the file is
%               generated with the same name as the border file, and
%               automatically deleted.
%   con_map :   Specifies which of the contrast maps to sample from
%               (default is all of the ones in gifti file)
%   b_name  :   If provided, only one border (specified by b_name) is used
%               for cross-section. By default, each border in border file
%               is used.
%
% OUTPUT:
%   Y       :   sampled data along the borders vertices
%
% To get the coordinates for the sampled points, sample also from the
% respective coord file.
% _________________________________________________________________
% Adapted from caret_crosssection (Joern Diedrichsen 2013)
% Last change: GA - 2019.05.15

% Defaults
outname = [];
con_map = [];
b_name  = [];
vararginoptions(varargin, {'outname', 'con_map', 'b_name'});

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

% Perform cross-section
if ~isempty(con_map)
    Y = D.cdata(B.cdata==1, con_map);
else
    Y = D.cdata(B.cdata==1, :);
end

% Delete border metric file
if saveout==0
    delete(outname);
end