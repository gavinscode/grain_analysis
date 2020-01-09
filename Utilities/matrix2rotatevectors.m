function [ rotationMatrix ] = matrix2rotatevectors( vectorA, vectorB )
% vector rotation
% Will calculate a matrix that rotates vector B to align with vector A

    % Correction orientation of vectors
    if size(vectorA,2) == 1, vectorA = vectorA'; end
    
    if size(vectorB,2) == 1, vectorB = vectorB'; end
    
    % Convert to unit vectors.
    vectorA = vectorA/norm(vectorA);
    
    vectorB = vectorB/norm(vectorB);
        
    isVectorsEqual = sum(vectorA == vectorB) == 3;
    isVectorsOpposite = sum(vectorA == -1*vectorB) == 3;
    
    if ~isVectorsEqual && ~isVectorsOpposite
        
        v = cross(vectorA, vectorB);
        
        rotationMatrix = eye(3) + calcSSC(v) + ...
            calcSSC(v)^2 * ( 1-dot(vectorA, vectorB) ) / ...
            ( norm( cross(vectorA, vectorB) )^2 );
        
    elseif isVectorsOpposite
            
        rotationMatrix = -1*eye(3);
        
        elseif isVectorsEqual
            
        rotationMatrix = eye(3);
        
    end
end

function ssc = calcSSC(v)
% Calculates skew symetric cross product of vector

    ssc = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];

end

