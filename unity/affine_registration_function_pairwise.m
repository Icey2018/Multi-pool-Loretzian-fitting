function [e]=affine_registration_function_pairwise(par,scale,Imoving,Ifixed,mtype,ttype,factor,interpMode)
% This function affine_registration_image, uses affine transfomation of the
% 3D input volume and calculates the registration error after transformation.
%
% I=affine_registration_image(parameters,scale,I1,I2,type);
%
% input,
%   parameters (in 2D) : Rigid vector of length 3 -> [translateX translateY rotate]
%                        or Affine vector of length 7 -> [translateX translateY  
%                                           rotate resizeX resizeY shearXY shearYX]
%
%   parameters (in 3D) : Rigid vector of length 6 : [translateX translateY translateZ
%                                           rotateX rotateY rotateZ]
%                       or Affine vector of length 15 : [translateX translateY translateZ,
%                             rotateX rotateY rotateZ resizeX resizeY resizeZ, 
%                             shearXY, shearXZ, shearYX, shearYZ, shearZX, shearZY]
%   
%   scale: Vector with Scaling of the input parameters with the same lenght
%               as the parameter vector.
%   I1: The 2D/3D image which is affine transformed
%   I2: The second 2D/3D image which is used to calculate the
%       registration error
%   mtype: Metric type: s: sum of squared differences.
%
% outputs,
%   I: An volume image with the registration error between I1 and I2
%
% example,
%
% Function is written by D.Kroon University of Twente (July 2008)
x = par .* scale;

% get the affine transformation matrix
M = transformation_matrix(ttype, x);

%  apply transformation of the moving img (3 stands for cubic interpolation)
I3 = affine_transform_2d_double(double(Imoving),double(M),interpMode); 

% metric computation
switch mtype
    case 'ssd' %squared differences
        e = sum((I3(:) - Ifixed(:)).^2) / numel(I3);
    case 'cc' % normalized intensity cross-correlation
        I3_mean = mean(I3(:));
        Ifixed_mean = mean(Ifixed(:));
        nom = sum( (I3(:) - I3_mean) .* (Ifixed(:) - Ifixed_mean) );
        denom = sqrt(sum((I3(:) - I3_mean).^2) .* sum((Ifixed(:) - Ifixed_mean).^2)) ;
        e = factor * nom / denom;
    case 'gcc' % normalized gradient cross-correlation
        [I3_fx,I3_fy] = imgradientxy(I3);
        [Ifixed_fx,Ifixed_fy] = imgradientxy(Ifixed);
        tmp = I3_fx .* Ifixed_fx + I3_fy .* Ifixed_fy;
        nom = sum(tmp(:));
        tmp = (I3_fx).^2 + (I3_fy).^2;
        tmp2 = (Ifixed_fx).^2 + (Ifixed_fy).^2;
        denom = sqrt(sum(tmp(:)) .* sum(tmp2(:))) ;
        e = factor * nom / denom;
    case 'mi'
        range = getrangefromclass(Ifixed);
        % make a joint image histogram
        bins = 32;
        [hist12, hist1, hist2] = mutual_histogram_double(double(Ifixed), double(I3), double(range(1)),double(range(2)),double(bins));
        
%         hist12 = imgaussian(hist12,1); hist1 = imgaussian(hist1,1); hist2 = imgaussian(hist2,1);
        
        %calculate probablities
        p1 = double(hist1); p2 = double(hist2); p12 = double(hist12);
       p1log=p1 .* log(p1+eps); p2log=p2 .* log(p2+eps); p12log=p12.* log(p12+eps);
         
       %calculate amount of information
       HA = -sum(p1log);  HB = -sum(p2log); HAB = -sum(p12log(:));
       
       %normalized mutual information
       if(HAB == 0), HAB = eps; end
       M = HA + HB - HAB;
       e = factor * M;         
    otherwise
        error('Unknown metric type');
end

