function [meanPofile, xCoordinate, stdProfile, countProfile, coordBins] = calculateProfile(values, coordinates, nBins, offset)

    % If offset is not empty, will offset to zero
    if isempty(offset)
        coordinates = coordinates - min(coordinates);
    else
        coordinates = coordinates - offset; 
    end

    if length(nBins) == 1
        % Add 1 as these are bin bounds
        binStep = floor((max(coordinates)-min(coordinates))/(nBins+1));

        coordBins = min(coordinates):binStep:max(coordinates);

        % Add just final point as ti will always be a bit more
        coordBins(end) = max(coordinates);
        
        xCoordinate = coordBins(1:end-1) + binStep/2;
    else
        coordBins = nBins;
        
        xCoordinate = coordBins(1:end-1) + (coordBins(2)-coordBins(1))/2;
    end
    
    meanPofile = zeros(length(coordBins)-1,1)*NaN;
    
    stdProfile = zeros(length(coordBins)-1,1)*NaN;
    
    countProfile = zeros(length(coordBins)-1,1)*NaN;
    
    % Take average and std of values in each bin
    for iBin = 1:(length(coordBins)-1)
        inds = find(coordinates > coordBins(iBin) & coordinates <= coordBins(iBin+1));
        
        if ~isempty(inds)        
            meanPofile(iBin) = mean(values(inds));
            
            stdProfile(iBin) = std(values(inds));
            
            countProfile(iBin) = length(inds);
        end
    end
end

