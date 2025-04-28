% function [regstered, tformSimilarity, Rfixed] = coRegister(fixS0, moveS0,transType)
function [movingRegistered] = coRegister(fixed, moving,transType)
%transType: 'rigid'/ 'affine'/ 'similarity'

[optimizer,metric] = imregconfig("multimodal");
% movingRegistered = imregister(moving,fixed,transType,optimizer,metric);
% 
% optimizer.InitialRadius = optimizer.InitialRadius/3.5;
% moving = imregister(moving,fixed,transType,optimizer,metric);
% % 
% optimizer.MaximumIterations = 100;
% moving = imregister(moving,fixed,transType,optimizer,metric);
% % 
% tformSimilarity = imregtform(moving,fixed,transType,optimizer,metric);
% % Rfixed = imref2d(size(fixed));
% % movingRegisteredSimilarity = imwarp(moving,tformSimilarity,OutputView=Rfixed);

% tformSimilarity =
% imregtform(moving,fixed,transType,optimizer,metric);%xuxi comment
% movingRegistered = imregister(moving,fixed,transType,optimizer,metric, ...
%     InitialTransformation=tformSimilarity); % original

movingRegistered = imregister(moving,fixed,transType,optimizer,metric); %xuxi use

% optimizer.InitialRadius = optimizer.InitialRadius/5;
% optimizer.InitialRadius = optimizer.InitialRadius/4.5;
% movingRegistered = imregister(moving,fixed,transType,optimizer,metric);

% optimizer.MaximumIterations = 100;
% movingRegistered = imregister(moving,fixed,transType,optimizer,metric);
end