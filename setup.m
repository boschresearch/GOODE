% Initialize Matlab for this project
% -------------------------------------------------------------------
% Copyright (c) 2019 Robert Bosch GmbH
% All rights reserved.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
% 
% Authors: 
%    David John (david.john@de.bosch.com)
%    Michael Schober (michael.schober5@de.bosch.com)

% Setup paths

% External GIT repos
addpath(genpath(fullfile(pwd, '..', '..', 'external', 'matlab2tikz', 'src')));

% Subdirectories
addpath(genpath(fullfile(pwd, 'external', 'bvpTestSet', 'matlabsrc')));
addpath(fullfile(pwd, 'external', 'tom'));

addpath(fullfile(pwd, 'experiments'));
addpath(genpath(fullfile(pwd, 'GPs')));
