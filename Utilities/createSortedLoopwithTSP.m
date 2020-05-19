function [coordinateListSorted, indexSorted] = createSortedLoopwithTSP(coordinateList)

    %3D part removed, but easy to add back in distance calcualtion.

    %so this copies a matlab optimisation toolbox examples of solving the travelling salesman problem.
    %talk about an indirect path to a solution!
    
    %note that this is completed based on 3D straight line distances, not geodesic paths...

    %reducer moved to calling code
%     toKeep = [];
%     toTest = 1:length(coordinateList);
%     while ~isempty(toTest);
%         dists = sqrt((coordinateList(toTest(1),1)-coordinateList(toTest,1)).^2 + (coordinateList(toTest(1),2)-coordinateList(toTest,2)).^2);
% 
%         toKeep = [toKeep toTest(1)];
%         toTest(dists < 2) = [];
%     end
%     coordinateList = coordinateList(toKeep, :);
    
    nStops = size(coordinateList,1); % you can use any number, but the problem size scales as N^2
    
%     if nStops > 100
%         warning(' a lot of points to solve');
%     end
    
    idxs = nchoosek(1:nStops,2);

    dist = sqrt((coordinateList(idxs(:,1),1) - coordinateList(idxs(:,2),1)).^2 + ...
                 (coordinateList(idxs(:,1),2) - coordinateList(idxs(:,2),2)).^2);
    lendist = length(dist);

    Aeq = spones(1:length(idxs)); % Adds up the number of trips
    beq = nStops;

    Aeq = [Aeq;spalloc(nStops,length(idxs),nStops*(nStops-1))]; % allocate a sparse matrix
    for ii = 1:nStops
        whichIdxs = (idxs == ii); % find the trips that include stop ii
        whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
        Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
    end
    beq = [beq; 2*ones(nStops,1)];

    intcon = 1:lendist;
    lb = zeros(lendist,1);
    ub = ones(lendist,1);

    opts = optimoptions('intlinprog','Display','off');
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

    tours = detectSubtours(x_tsp,idxs);
    numtours = length(tours); % number of subtours

    A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
    b = [];
    its = 0;
    while numtours > 1 % repeat until there is just one subtour
        % Add the subtour constraints
        b = [b;zeros(numtours,1)]; % allocate b
        A = [A;spalloc(numtours,lendist,nStops)]; % a guess at how many nonzeros to allocate
        for ii = 1:numtours
            rowIdx = size(A,1)+1; % Counter for indexing
            subTourIdx = tours{ii}; % Extract the current subtour
    %         The next lines find all of the variables associated with the
    %         particular subtour, then add an inequality constraint to prohibit
    %         that subtour and all subtours that use those stops.
            variations = nchoosek(1:length(subTourIdx),2);
            for jj = 1:length(variations)
                whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                           (sum(idxs==subTourIdx(variations(jj,2)),2));
                A(rowIdx,whichVar) = 1;
            end
            b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
        end

        % Try to optimize again
        [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);

        % Visualize result
%         lh = updateSalesmanPlot(lh,x_tsp,idxs,stopsLon,stopsLat);

        % How many subtours this time?
        tours = detectSubtours(x_tsp,idxs);
        numtours = length(tours); % number of subtours
        its = its + 1;
        
        if its > 3
            warning('line sorter iterating too long, probably multi region - its %i', its);
        end
    end
    
    segments = find(x_tsp); % Indices to trips in solution
    coordinateListSorted = zeros(size(coordinateList));
    indexSorted = zeros(size(coordinateList,1),1);
    
    coordinateListSorted(1,:) = coordinateList(idxs(segments(1),1),:);
    indexSorted(1) = idxs(segments(1),1);
    
    starts = idxs(segments(:),1);
    stops = idxs(segments(:),2);
    
    prev = stops(1);
    lastIsStop = 1;
    linesNotUsed = ones(length(segments), 1);
    linesNotUsed(1) = 0;
    
%     [starts stops]
    
    %fml, sorting this was not intuitive...
    
%   figure; hold on;
    for i = 1:length(segments)-1
        if i >= length(segments)
            break;
        end
        if lastIsStop
            ind = find(starts == prev & linesNotUsed);
            if isempty(ind)
                ind = find(stops == prev & linesNotUsed);
                lastIsStop = 0;
                if isempty(ind);
                    warning('line sorter may have problems');
                    break;
                end
                prev = starts(ind(1));
            else
                prev = stops(ind(1));
            end
        else
            ind = find(stops == prev & linesNotUsed);
            if isempty(ind)
                ind = find(starts == prev & linesNotUsed);
                lastIsStop = 1;
                if isempty(ind);
                    warning('line sorter may have problems');
                    break;
                end
                prev = stops(ind(1));
            else
                prev = starts(ind(1));
            end
        end
        ind = ind(1);
%         [ind prev lastIsStop]

        
        linesNotUsed(ind) = 0;

%       line([coordinateList(starts(i),1) coordinateList(stops(i),1)], [coordinateList(starts(i),2) coordinateList(stops(i),2)]);
        if lastIsStop
            coordinateListSorted(i+1,:) =  coordinateList(starts(ind),:);
            indexSorted(i+1) = starts(ind);
        else
            coordinateListSorted(i+1,:) =  coordinateList(stops(ind),:);
            indexSorted(i+1) = stops(ind);
        end
    end
        
    if size(unique(coordinateListSorted,'rows'),1) ~= length(segments)
        error('problem with line sorter');
    end
%         figure;
%         plot3(coordinateListSorted(:,1), coordinateListSorted(:,2));
end