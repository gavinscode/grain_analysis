function [normal] = normaltosurface( pointCoords, surface, referenceCenter, flipAway, includeRange) 
    % Flip away will flip eventual normal away from reference centre,
        % other wise normal will be flipped towards reference centre.
    % Reference centre should be rough direction of inside of convex strucutre (if availible.)

    % As in bee code, however, radius adjust from levels sets and height req. removed.

    debugPlot = 0;
    
    if isempty(includeRange)
        includeRange = 50;   % A bit arbitary...
    end
    
    if debugPlot 
        figure;
    end

    % Get indices to use and original point
    surfIndexs = find( sqrt( (surface(:,1)-pointCoords(1)).^2 + ...
        (surface(:,2)-pointCoords(2)).^2 + (surface(:,3)-pointCoords(3)).^2) < includeRange);
        
    [~, orignalPointInd] = min( sqrt( (surface(surfIndexs,1)-pointCoords(1)).^2 + ...
        (surface(surfIndexs,2)-pointCoords(2)).^2 + (surface(surfIndexs,3)-pointCoords(3)).^2));
    
    if debugPlot
        subplot(1,2,1); hold on; axis equal
        plot3(surface(surfIndexs,1), surface(surfIndexs,2), surface(surfIndexs,3), 'r.');
    end
    
    % Rotate surface to axis of max variation.
    pointsCenter = mean(surface(surfIndexs,:));
    
    surfPoints = surface(surfIndexs,:)-pointsCenter;
    
    pcaVecs = pca(surfPoints);

    %vectorRotationFromAtoB([0 0 1]',pcaVecs(:,3)');
    forwardsRotation = matrix2rotatevectors([0 0 1]',pcaVecs(:,3)');
    
    backwardsRotation = matrix2rotatevectors(pcaVecs(:,3)',[0 0 1]');  
        
    % Equivlent to test below.
    if ~isempty(referenceCenter)
        referenceCenter = referenceCenter-pointsCenter; 

        referenceCenter = referenceCenter*forwardsRotation;
    end
    
    % Test if ref center points down, if so rotate up.
    %%% Removed section as convex no longer forced
%     referenceCenterOriginal = referenceCenter;
%     
%     referenceCenter = referenceCenterOriginal-pointsCenter; 
%     
%     referenceCenter = referenceCenter*forwardsRotation;
%     
%     if referenceCenter(3) < 0
%         
%         forwardsRotation = matrix2rotatevectors([0 0 -1]',pcaVecs(:,3)');
%         
%         backwardsRotation = matrix2rotatevectors(pcaVecs(:,3)',[0 0 -1]');
%         
%         referenceCenter = referenceCenterOriginal-pointsCenter;
%         
%         referenceCenter = referenceCenter*forwardsRotation;
%         
%     else
%         
%         backwardsRotation = matrix2rotatevectors(pcaVecs(:,3)',[0 0 1]');  
%     end

    newPoints = surfPoints*forwardsRotation;    
            
    % Shift base to zero as cf tool can't translate horizontally.
    [~, minimaInd] = min(newPoints(:,3));
    
    originOffset = newPoints(minimaInd,:);
    
    newPoints = newPoints - originOffset; 

    if debugPlot
         subplot(1,2,2); hold on; axis equal
         plot3(newPoints(:,1), newPoints(:,2), newPoints(:,3), 'r.','markersize',20);
         if ~isempty(referenceCenter)
            line([0 referenceCenter(1)], [0 referenceCenter(2)], [0 referenceCenter(3)], 'color', 'r');
         end
    end
    
    opts = fitoptions( 'Method', 'LinearLeastSquares' );

    % Can force convex, but not always true.
    %opts.Lower = [-Inf -Inf -Inf 0 0 0]; 
    %opts.Upper = [Inf Inf Inf Inf Inf Inf];

    fittedPoly = fit(newPoints(:,1:2), newPoints(:,3),'poly22', opts);
     
    if debugPlot
        zlim([-100 100]); xlim([-500 500]); ylim([-500 500]);
        
        plot(fittedPoly);
    end
    
    % Find normal at point that original was shifted to
    pointToTest = newPoints(orignalPointInd,:);
    
    [dX, dY] = differentiate(fittedPoly, pointToTest(1), pointToTest(2));

    % Correct Nan values to zero.
    if isnan(dX); dX = 0; end
    
    if isnan(dY); dY = 0; end

    % Compose normal from derivatives.
    N = [dX, dY, -1]; 
    
    % Normalize normal and reference centre vectors.
    N = N/sqrt(N(1)^2+N(2)^2+N(3)^2);
    
    if ~isempty(referenceCenter)
        referenceCenter = referenceCenter/sqrt(referenceCenter(1)^2+referenceCenter(2)^2+referenceCenter(3)^2);

        % Flip normal to align with reference centre
        if ~flipAway

            if sqrt((referenceCenter(1)-N(1))^2+(referenceCenter(2)-N(2))^2+(referenceCenter(3)-N(3))^2) > sqrt(2)

                N = N*-1;
            end

        else

            if sqrt((referenceCenter(1)-N(1))^2+(referenceCenter(2)-N(2))^2+(referenceCenter(3)-N(3))^2) < sqrt(2)

                N = N*-1;
            end

        end
    end
    
    if debugPlot
        line([pointToTest(1) pointToTest(1)+N(1)*100],[pointToTest(2) pointToTest(2)+N(2)*100] , [0 N(3)*100]);
    end
    
    normal = N*backwardsRotation;

    if debugPlot
         subplot(1,2,1);
         line([pointCoords(1) pointCoords(1)+normal(1)*50],[pointCoords(2) pointCoords(2)+normal(2)*50] , [pointCoords(3) pointCoords(3)+normal(3)*50]);
    end
end
