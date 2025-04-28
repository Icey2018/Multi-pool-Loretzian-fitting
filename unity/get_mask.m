function [Mask_use] = get_mask(Mask)
% make multiple Mask into one
% input: mask:[x,y,num]
% output: Mask_use:[x,y]
Mask_use = Mask(:,:,1);
for s = 2:size(Mask,3)
    Mask_use = Mask_use-Mask(:,:,s);
end

end