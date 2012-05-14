%This script is used to load the data necessary for running the genesis
%potential index algorithm

tic
if ~exist('dataLoaded', 'var') || ~dataLoaded
    dataPath = '/project/expeditions/lem/data/pressureLevelData_1979-present.nc';
    
    lat = ncread(dataPath, 'lat');
    lon = ncread(dataPath, 'lon');
    
    time = ncread(dataPath, 'time');
    
    levels = ncread(dataPath, 'lev');
    %change levels from Pa to mb
    levels = levels*.01;
    
    relHumidity = ncread(dataPath, 'var157');
    %only need relative humidity from 700mb
    relHumidity = relHumidity(:, : , levels(:) == 700, :);
    relHumidity = permute(relHumidity, [2 1 4 3]);
    
    relVorticity = ncread(dataPath, 'var138');
    %only need vorticity at 850 mb
    relVorticity = relVorticity(:, :, levels(:) == 850, :);
    relVorticity = permute(relVorticity, [2 1 4 3]);

    earthRotation = 2*pi/24/60/60; %rad/sec
    %earth vorticity assumes that degrees are always positive.  i.e. they
    %range from 90 degrees north to 90 degrees south
    earthVorticity = 2 * earthRotation * sin(degtorad(abs(lat)));
    absVorticity = zeros(size(lat, 1), size(lon, 1), size(time, 1));
    for i = 1:size(lon, 1)
        for j = 1:size(time, 1)
            absVorticity(:, i, j) = earthVorticity + relVorticity(:, i, j);
        end
    end

    windSpeeds = ncread(dataPath, 'var131');
    %calculate the vertical wind shear between 850 and 200 mb
    vWindShear = windSpeeds(:, :, levels(:) == 850, :) - windSpeeds(:, :, levels(:) == 200, :);
    vWindShear = permute(vWindShear, [2 1 4 3]);
    load PIMaps.mat;
    maxWindSpeeds = zeros(size(lat, 1), size(lon, 1), size(time, 1));
    for i = 1:size(time, 1)
        maxWindSpeeds(:, :, i) = PIData{i, 1};
    end
    
    dataLoaded = true;
end

GPIData = cell(size(time, 1), 1);
parfor currTime = 1:size(time, 1)
    gpiMat = zeros(256, 512);
    for i = 1:size(lat, 1)
        gpiRow = zeros(1, size(lon, 1));
        for j = 1:size(lon, 1)
            if maxWindSpeeds(i, j, currTime) == 0
                continue;
            end
            gpiRow(1, j) = gpi(absVorticity(i, j, currTime),relHumidity(i, j, currTime), maxWindSpeeds(i, j ,currTime),vWindShear(i, j, currTime));
        end
        gpiMat(i, :) = gpiRow;
    end
    currTime
    GPIData{currTime} = gpiMat;
end
    
time = toc;
save('GPIData.mat', 'GPIData', 'time');


