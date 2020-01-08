function [ arrayOut ] = loadtiffstack( address, loadSmall)
% Copy tiff stack saved by amira into volume. 
    %colour image will be converted to greyscale, but unique values for
    %colours are not guaranteed
% address - is folder to load
% loadSmall - flag, if set will arrayOut will be in uint8 format, 
              %otherwise, will be double
                   
    currentDirectory = pwd;          
              
    % Find which files to load.
    cd(address); dirtoryContents = dir; includedFiles = [];
    
    for i = 1:length(dirtoryContents)
        
        if strfind(dirtoryContents(i).name, 'tif')   
            
            includedFiles = [ includedFiles i ];
            
        end
        
    end
    
    
    
    % Load each file into an array.
    counter = 1; nImages = length(includedFiles);
    
    for  i = 1:nImages
        
        % Initalize array given dimensions of first image.
        if i == 1
            
            [tempImage] = imread(dirtoryContents(includedFiles(i)).name, 'TIFF');
            
            if loadSmall
                
                arrayOut = zeros(size(tempImage,1), size(tempImage,2), nImages, 'uint8');
                
                arrayOut(:,:,counter) = uint8(tempImage);
                
            else
                
                arrayOut = zeros(size(tempImage,1), size(tempImage,2), nImages);
                
                arrayOut(:,:,counter) = tempImage;
                
            end
            
        else
            
            tempImage = imread(dirtoryContents(includedFiles(i)).name, 'TIFF');
            
            if loadSmall
                
                arrayOut(:,:,counter) = uint8(tempImage);
                
            else
                
                arrayOut(:,:,counter) = tempImage;
                
            end
            
        end
        
        counter = counter + 1;
    
    end
    
    cd(currentDirectory)
    
end

