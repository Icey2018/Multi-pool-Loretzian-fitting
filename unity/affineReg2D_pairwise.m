function [Iregistered, x] = affineReg2D_pairwise(Im,If, mtype, ttype, x, scale, factor, interpMode, plotMetric)
% Example of 2D affine registration
%   Robert Mart? (robert.marti@udg.edu)
%   Based on the files from  D.Kroon University of Twente 

x = x ./ scale;

% optimizer  
opts = optimset('Display', 'off', ...
                'MaxIter', 3000, ...
                'TolFun', 1.0e-10, ...
                'TolX', 1.0e-10, ... 
                'MaxFunEvals', 10000*length(x));
                %'PlotFcns',@optimplotfval,...
if (plotMetric)
    opts = optimset(opts,'PlotFcns',@optimplotfval);                
end

[x] = fminsearch(@(x)affine_registration_function_pairwise(x,scale,Im,If,mtype,ttype,factor,interpMode),...
    x,opts);
% fminunc

x = x .* scale;

% get the transformation matrix
M = transformation_matrix(ttype, x);

 % Transform the image (3 stands for cubic interpolation)
Iregistered = affine_transform_2d_double(double(Im),double(M),interpMode); 


end

