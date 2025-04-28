function Iregistered = image_registration_pairwise (Imoving,Ifixed, mtype, ttype,depth)
%Input:
%Imiving: size: Nx,Ny
%Ifixed: sixe: Nx,Ny
% mtype: metric type
% 'ssd': ssd | 'cc': cross-correlation | 'gcc': gradient correlation |
% 'mi': mutual inforamtion
% ttpye: rigid registration
% 'rigid': rigid | 'affine': affine (7-params) | 'affine6params': affine (6-params)
%depth: multiresolution levels (5 is enough)
%output: Iregistered image
tStart = tic;

% if need to display "metric vs iteration" plot
plotMetric = false;
interpMode = 0;%3;

factor = -1; % change to 1 if using images with inverted intensities like brain4


% build the pyramids in a decreasing, bottom-up order
% declare the pyramids
ImovingPyr = cell(depth);
IfixedPyr = cell(depth);
% temporary aliases
Imoving_ = Imoving;
Ifixed_ = Ifixed;
for i = 1:depth
    % add the images to their pyramids
    ImovingPyr{depth-i+1} = Imoving_;
    IfixedPyr{depth-i+1} = Ifixed_;
    % subsample
    Imoving_ = impyramid(Imoving_, 'reduce');
    Ifixed_ = impyramid(Ifixed_, 'reduce');
end

% initial variables values
switch ttype
    case 'rigid'
        params = [0 0 0];
        scale = [1 1 0.1];
    case 'affine6params'
        params = [0  0     1     0     0     1];
        scale  = [1  1 0.001 0.001 0.001 0.001];
    case 'affine'
        params = [0  0    0   1   1      0      0];
        scale  = [1  1 0.01 0.1 0.1 0.0001 0.0001];
end

% iteratively register images at different levels of the pyramid
for i = 1:depth
    % refer to affine_transform_2d_double.m for the different options
%     if i == depth
%         interpMode = 4;
%     else
%         interpMode = 2;
%     end
    
%     fprintf('Registering pyramids at level %i:\n', i);
    Imoving_ = ImovingPyr{i};
    Ifixed_ = IfixedPyr{i};
    [Iregistered, x] = affineReg2D_pairwise(Imoving_, Ifixed_, mtype, ttype, params, scale, factor, interpMode, plotMetric);
    x(1) = x(1) * 2;
    x(2) = x(2) * 2;
    params = x;
end

tElapsed = toc(tStart);


