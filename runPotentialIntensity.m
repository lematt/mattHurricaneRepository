

%This script is used to load the data necessary to run the potential intensity algorithm mpikerry.m.

tic

if ~exist('dataLoaded', 'var') || ~dataLoaded
    % temps is a lon x lat x pressure level x time, matrix
    temps = ncread('/project/expeditions/lem/data/temperature_eraInterim_1979-present.nc', 'var130');
    
    % temps is a lat x lon x time x pressure level,  matrix
    temps = permute(temps, [2 1 4 3]);
    
    % change temperature from kelvin to celsius 
    temps = temps - 273.15;
    
    % load list of times the data was taken at
    time = ncread('/project/expeditions/lem/data/temperature_eraInterim_1979-present.nc', 'time');
    
    % load latitude and longitude
    lat = ncread('/project/expeditions/lem/data/temperature_eraInterim_1979-present.nc', 'lat');
    lon = ncread('/project/expeditions/lem/data/temperature_eraInterim_1979-present.nc', 'lon');
    
    % load the list of pressure levels
    levels = ncread('/project/expeditions/lem/data/temperature_eraInterim_1979-present.nc', 'lev');
    [sortedLevels, sortIndices] = sort(levels, 'descend');
    sortedLevels = double(sortedLevels);
    
    %load the specific humidity data
    specificHumidity = ncread('/project/expeditions/lem/data/specificHumidity_eraInterim_1979-present.nc', 'var133');
    specificHumidity = permute(specificHumidity, [ 2 1 4 3 ]);
    
    % Verify that the levels are the same in specific humidity data
    splevels = ncread('/project/expeditions/lem/data/specificHumidity_eraInterim_1979-present.nc', 'lev');
    if ~all(splevels == levels)
        error('levels in specific humidity are not the same as in temperature')
    end
    
    %change levels from Pa to mb
    levels = levels * .01;
    
    % Reverse the order of the pressure levels in the spec. hum. and temp.
    temps = temps(:, :, :, sortIndices);
    specificHumidity = specificHumidity(:, :, :, sortIndices);
    
    % calculate mixing ratio from specific humidity.  w = q/(1-q) where w =
    % mixing ratio and q = specific humidity
    mixingRatio = specificHumidity ./ (1 - specificHumidity);

    sst = ncread('/project/expeditions/lem/data/sst_slp_eraInterim_1979-2010.nc', 'var34');
    sst = permute(sst, [2 1 3]);
    sst = sst - 273.15;
    
    centralPressure = ncread('/project/expeditions/lem/data/sst_slp_eraInterim_1979-2010.nc', 'var134');
    centralPressure = permute(centralPressure, [2 1 3]);
    %change from Pascal to millibars
    centralPressure = centralPressure * .01;
    
    dataLoaded = true;
end

%reverse the order of pressure levels, they need to be from highest to
%lowest for the potential intensity algorithm
% n = size(levels, 1);
% pressureLevels = zeros(n, 1);
% 
% t=temps(1, 1, :, 1);
% mr = mixingRatio(1, 1, :, 1);
% for i = 1:n
%     pressureLevels(i) = levels(n - i + 1);
%     t(i) = temps(120, 240, n - i + 1, 1);
%     mr(i) = mixingRatio(120, 240, (n - i + 1), 1);
% end

%right now we are only trying to run the algorithm for one point on the map
%and for one point in time.

%t = permute(temps(120, 240, 1, :), [ 4 1 2 3 ]);
%r = permute(mixingRatio(120, 240, 1, :), [ 4 1 2 3 ]);

%takes sst, central pressure, pressure levels, temperature, and mixing
%ratio.  The last three argumenst are vectors.
% [pmin,vmax,capea,ifl] = mpikerry(sst(120, 240, 1), centralPressure(120, 240, 1), pressureLevels, t, mr)
%[pmin,vmax,capea,ifl] = mpikerry(sst(128, 256, 1), centralPressure(60, 70, 1), sortedLevels, t, r);
%pmin = minimum central pressure in mb.  vmax = maximum velocity, if ifl =
%0 or 1 (0 = no convergence, 1 = ok, and if ifl = 2 then algorithm failed

vmaxMap = zeros(256, 512);
pminMap = zeros(256, 512);

PIData = cell(size(time, 1), 2);

for currentTime = 1:size(time)
    for i = 1:size(lat)
        for j = 1:size(lon)
            t = permute(temps(i, j, currentTime, :), [4 1 2 3]);
            r = permute(mixingRatio(i, j, currentTime, :), [ 4 1 2 3]);
            [pmin, vmax, capea, ifl] = mpikerry(sst(i, j, currentTime), centralPressure(i, j, currentTime), sortedLevels, t, r);
            vmaxMap(i, j) = vmax;
            pminMap(i, j) = pmin;
        end
    end
    PIData{currentTime, 1} = vmaxMap;
    PIData{currentTime, 2} = pminMap;
end
%imagesc(vmaxMap)
save('PIMaps.mat', 'PIData');

toc
