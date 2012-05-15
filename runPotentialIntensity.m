

%This script is used to load the data necessary to run the potential intensity algorithm mpikerry.m.
%information regarding the data used can be found here http://rda.ucar.edu/datasets/ds627.1/metadata/grib.html
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
    
    % load the list of pressure levels so that they are indexed from the
    % surface upwards
    levels = ncread('/project/expeditions/lem/data/temperature_eraInterim_1979-present.nc', 'lev');
    [sortedLevels, sortIndices] = sort(levels, 'descend');
    sortedLevels = double(sortedLevels);
    
    %load the specific humidity data (data is in kg/kg)
    specificHumidity = ncread('/project/expeditions/lem/data/specificHumidity_eraInterim_1979-present.nc', 'var133');
    specificHumidity = permute(specificHumidity, [ 2 1 4 3 ]); 
    
    % Verify that the levels are the same in specific humidity data
    splevels = ncread('/project/expeditions/lem/data/specificHumidity_eraInterim_1979-present.nc', 'lev');
    if ~all(splevels == levels)
        error('levels in specific humidity are not the same as in temperature')
    end
    
    %change levels from Pa to mb
    levels = levels * .01;
    
    % Reverse the order of the pressure levels in the specific humidity and temp.
    temps = temps(:, :, :, sortIndices);
    temps = permute(temps, [4, 1, 2, 3]); %lev x lat x lon x time
    specificHumidity = specificHumidity(:, :, :, sortIndices);
    
    % calculate mixing ratio from specific humidity.  w = q/(1-q) where w =
    % mixing ratio and q = specific humidity
    mixingRatio = specificHumidity ./ (1 - specificHumidity);
    mixingRatio = permute(mixingRatio, [4, 1, 2, 3]); %levxlatxlonxtime

    sst = ncread('/project/expeditions/lem/data/sst_slp_eraInterim_1979-2010.nc', 'var34');
    sst = permute(sst, [2 1 3]);
    sst = sst - 273.15; %change from kelvin to celsius
    
    centralPressure = ncread('/project/expeditions/lem/data/sst_slp_eraInterim_1979-2010.nc', 'var134');
    centralPressure = permute(centralPressure, [2 1 3]);
    %change from Pascal to millibars
    centralPressure = centralPressure * .01;
    
    dataLoaded = true;
end

vmaxMap = zeros(256, 512);
pminMap = zeros(256, 512);

PIData = cell(size(time, 1), 2);            
currTime = 1;
for currentTime = 1:size(time)
    for i = 1:size(lat)
        for j = 1:size(lon)
            if isnan(sst(i, j, currTime)) %coordinate is on land, no sst
                continue;
            end
            [pmin, vmax, ~,~] = mpikerry(sst(i, j, currTime), centralPressure(i, j, currTime), sortedLevels, temps(:, i, j, currTime), mixingRatio(:, i, j, currTime));
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
