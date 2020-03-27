% Compute within- and between-event residuals
% Jack Baker
% Created 6/2/2016
% Updated 9/30/2019 for spatial correlation study, to revise data filtering
% and provide an end-to-end code base for analysis


clear; close all; clc;

%% load raw data
load NGA_W2_corr_meta_data % needed data from the NGA-West2 flatfile


%% initialize matrices

% trim long periods 
tIndex = find(Periods<=10); 
Periods = Periods(tIndex);
Sa_RotD50 = Sa_RotD50(:,tIndex);

% get sizes of matrices
numRecs = length(eqid);
numT = length(Periods);

[eventEQID,eventIdx] = unique(eqid); % get indices of unique events
eventEQID(1) = []; % throw out the first case, which is eqid = -999
eventIdx(1) = []; % throw out the first case, which is eqid = -999
for i = 1:length(eventEQID) % for each earthquake
   idx = find(eqid == eventEQID(i),1); % get index of first recording from this earthquake
   eventMag(i,1) = magnitude(idx); % get event magnitude of this earthquake
end

% initialize residual matrices
resid_RotD50Within  = nan*ones(numRecs, numT);
resid_RotD50Total   = nan*ones(numRecs, numT);

resid_RotD50BetweenLong  = nan*ones(max(eqid), numT); % leave a row for each eqid

% get max usable period for each ground motion
maxUsableT = 1./lowest_usable_freq;
% screening flag based on Vs30 (to match Heresi and Miranda)
allowableVs30 = (soil_Vs30 >= 180 &  soil_Vs30 <= 760);
% screening flag based on distance (to omit distant low-amplitude motions and missing-value cases)
allowableR = (closest_D < 300 & closest_D > 0);
allowableM = (magnitude >=4); % omit small-magnitude events



%% get predictions and Sa residuals

% get fw/hw terms for CY 2014 GMPE
FwHw = zeros(numRecs,1);
for i = 1:numRecs
    if strcmp(Fhw{i},'hw')
        FwHw(i) = 1;
    end
end

 
for i = 1:numT % period index
    i
    allowableFilter = (maxUsableT > Periods(i) & Sa_RotD50(:,i)>0); % omit missing records spectra that are outside usable period range or are missing
%     idxTemp = find(allowableFilter & chiouYoungsUsed); % find records used by Chiou and Youngs, and with usable data at this period
    idxTemp = find(allowableFilter & allowableVs30 & allowableR); % find records used by Chiou and Youngs, and with usable data at this period
    [~, eqIdx] = sort(eqid(idxTemp)); % now adjust index so that it puts the data in order by eqid (needed for mixed effects function)
    idx = idxTemp(eqIdx);

    numUsableRecs(i) = length(idx);
    
    %%%%%%%%%%%%%%%%%%%%%% RotD50 %%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:length(idx)
        [median_GMPE(k,1), sigma_GMPE(k,1), phi_GMPE(k,1), tau_GMPE(k,1)] = CY_2014_nga_mod(magnitude(idx(k)), Periods(i), closest_D(idx(k)), ...
                                                    Rjb(idx(k)), Rx(idx(k)), Z_tor(idx(k)), dip(idx(k)), rakeAngle(idx(k)), ...
                                                    Z1_CVMH(idx(k)), soil_Vs30(idx(k)), FwHw(idx(k)), 1, Region_BSSA(idx(k)));
    end
    
    % build data structure
    input.imObservations = Sa_RotD50(idx,i);
    input.imMedian   = median_GMPE;
    input.totalSigma = sigma_GMPE;
    input.interSigma = tau_GMPE;
    input.intraSigma = phi_GMPE;
    input.eventIds   = eqid(idx);
    Out=GetInterIntraEventResiduals(input);

    
    % store results
    resid_RotD50Total(idx,i)  = (log(Sa_RotD50(idx,i)) - log(median_GMPE))./sigma_GMPE;
    resid_RotD50Within(idx,i)  = Out.intraEventResidualsNormalized;
    medianPred(idx,i) = median_GMPE;
    
    idxNgms = find(Out.eventData.eventNumGms <= 1);
    Out.eventData.interEventResidualNormalized(idxNgms) = nan;

    resid_RotD50BetweenLong(Out.eventData.eventId ,i) = Out.eventData.interEventResidualNormalized; % alternate set of event residuals (for debugging)
    resid_BetweenNumRecs(Out.eventData.eventId ,i) = Out.eventData.eventNumGms;

    clear median_GMPE sigma_GMPE phi_GMPE tau_GMPE % clean up variables before moving to the next period
end



%% get event magnitudes for the between event residuals
eventMagLong = nan*ones(max(eqid),1);
for i = 1:length(eventEQID)
    eventMagLong(eventEQID(i)) = eventMag(i);
end


%% save data

save('mixedEffectsResids.mat', 'resid_RotD50Within',  'resid_RotD50BetweenLong',  'resid_RotD50Total',  'medianPred', ...
                                'resid_BetweenNumRecs', 'Periods', 'eventMagLong', 'eqid', 'magnitude', 'soil_Vs30', 'Rjb', 'numUsableRecs')

%% plots
% number of records per period
figure
semilogx([Periods ], [numUsableRecs ], '-b')
% axis([0.01 11 0 15000])
xlabel('Period (s)')
ylabel('Number of ground motions')
FormatFigure



