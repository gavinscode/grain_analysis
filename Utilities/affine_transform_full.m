function Iout=affine_transform_full(Iin,M,mode)
% Function affine_transform, is a wrapper of the (mex) functions
% affine_transform_2d_double and affine_transform_3d mex functions
%
% Iout = affine_transform(Iin,M,mode)
%
% inputs,
%   Iin :  Input image.
%   M : Transformation matrix
%   mode: If 0: linear interpolation and outside pixels set to nearest pixel
%            1: linear interpolation and outside pixels set to zero
%            2: cubic interpolation and outside pixels set to nearest pixel
%            3: cubic interpolation and outside pixels set to zero
%            (4): nearest interpolation and outside pixels set to nearest pixel - not added for 3D
%            5: nearest interpolation and outside pixels set to zero
% outputs,
%   Iout: Output image

% Function is written by D.Kroon University of Twente (September 2008)

% Modified by Gavin Taylor to prevent 2D execution. (Jan 2020) 
% Nearest neihbour interpolation was enabled in following c files:
% affine_transform_3d_double & affine_transform_3d_single & image_interpolation
% these are placed in original directory via symbolic links for compilation
% 
% To do: Add support for nearest pixel on border (4) by adding:
%   interpolate_3d_float_nearest & interpolate_3d_nearest to image_interpolation
% Also add/check support for 2D

if(~exist('mode','var')), mode=0; end
if(size(Iin,3)<4)
    if(~isa(Iin,'double')), Iin=im2double(Iin); end
    error('Not implemented in 2d');
	Iout=affine_transform_2d_double(double(Iin),double(M),double(mode));
else
    if(isa(Iin,'double'))
		Iout=affine_transform_3d_double(double(Iin),double(M),double(mode));
    else
		Iout=affine_transform_3d_single(single(Iin),single(M),single(mode));
    end
end
