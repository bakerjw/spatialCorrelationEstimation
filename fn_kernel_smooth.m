function [xSort,yVals] = fn_kernel_smooth(xObs, yObs)
% simple windowed smoothing of data, to help with plotting of trends
%
% Created by Jack Baker, 10/7/2019

nAvg = 25; % how many nearest neighbors to average

% sort the arrays into ascending order
[xSort,idx] = sort(xObs);
ySort = yObs(idx);

yVals = movmean(ySort,nAvg);

end

