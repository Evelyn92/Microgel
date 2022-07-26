function screen2tif(filename)
%SCREEN2JPEG Generate a JPEG file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
%
%   SCREEN2JPEG('filename') saves the current figure to the
%   JPEG file "filename".
%
%    Sean P. McCarthy
%    Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved
if nargin < 1
     error('Not enough input arguments!')
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
print -dtiff filename.tiff -r600