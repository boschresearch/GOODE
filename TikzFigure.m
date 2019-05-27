function TikzFigure(filename, fig_path)
%TikzFigure Helper function to plot with matlab2tikz

% How to use matlab tikz plot
% addpath(genpath(fullfile(pwd, '..', 'matlab2tikz', 'src')));
% t = linspace(0,1,10);
% y = sin(t);
% figure()
% plot(t,y)
% TikzFigure('testfigure.tex')
%
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

if nargin < 2
    fig_path = 'doc/icml2019/figures/';
end

filename = [fig_path filename];

x0=300;
y0=130;
width=550;
height=400;
set(gcf,'position',[x0,y0,width,height])
try
    cleanfigure()
    matlab2tikz(filename, 'height', '\figureheight', 'width', '\figurewidth');
catch e
    fprintf(2,'There was an error! The message was: %s\n',e.message);
end

end