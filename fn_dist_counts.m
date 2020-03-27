function [cumPairs] = fn_dist_counts (lats, longs, hFine)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jack Baker 2/1/2019
% last updated 10/2/2019
%
% Compute number of station pairs within each distance, for each earthquake
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%   lats = list of site lats
%   longs = list of site longs
%   options.maxR = maximum distance to which the variogram is computed
%
% Output Variables
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Compute site-to-site distances
nsites=length(lats);
        
distance = zeros(nsites);
for i=1:nsites
    for j=1:nsites
        distance(i,j) = pos2dist(lats(i),longs(i),lats(j),longs(j),1); 
    end
end

[cumPairs] = histcounts(distance(:),[0 hFine],'Normalization','cumcount');


end