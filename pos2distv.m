function dist = pos2distv(XN, XM)

% Vectorized version by Lukas Bodenmann, 04/13/2022
% Supports currently only method 1: plane approximation,
% which was also used in the original study.
% XN and XM are matrices with latitudes in the first column and longitudes in the
% second column.
% The result is a matrix with distances in km between each (lon,lat) pair
% in XN and each (lon,lat) pair in XM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate distance between two points on earth's surface
% given by their latitude-longitude pair.
% Input lag1,lon1,lag2,lon2 are in degrees, without 'NSWE' indicators.
% Input method is 1 or 2. Default is 1.
% Method 1 uses plane approximation,
% only for points within several tens of kilometers (angles in rads):
% d =
% sqrt(R_equator^2*(lag1-lag2)^2 + R_polar^2*(lon1-lon2)^2*cos((lag1+lag2)/2)^2)
% Method 2 calculates sphereic geodesic distance for points farther apart,
% but ignores flattening of the earth:
% d =
% R_aver * acos(cos(lag1)cos(lag2)cos(lon1-lon2)+sin(lag1)sin(lag2))
% Output dist is in km.
% Returns -99999 if input argument(s) is/are incorrect.
%
% Flora Sun, University of Toronto, Jun 12, 2004.
%
% Downloaded from Matlab File Exchange, 
% https://www.mathworks.com/matlabcentral/fileexchange/5256-pos2dist
%
% Copyright 2004, Flora Sun
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

lats1 = XN(:,1);
longs1 = XN(:,2);
longs1(longs1<0) = longs1(longs1<0) + 360;
XM = XM';
lats2 = XM(1,:);
longs2 = XM(2,:);
longs2(longs2<0) = longs2(longs2<0) + 360;

km_per_deg_la = 111.3237;
km_per_deg_lo = 111.1350;
km_la = km_per_deg_la * (lats1 - lats2);
dif_lo = abs(longs1 - longs2);
dif_lo(dif_lo > 180) = dif_lo(dif_lo > 180) - 180;
km_lo = km_per_deg_lo * dif_lo .* cos((lats1+lats2).*pi./360);
dist = sqrt(km_la.^2 + km_lo.^2);