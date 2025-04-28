function [z_B0correction, B0map, delta_freq_map]= B0correction(Para, images,Mask)
% B0map: x,y,slice,1,b1
% delta_freq_map: x,y,slice,1,b1
% z_B0correction: x,y,slice,freq, b1
%  interpolate step from 0.2 ppm to 0.1 ppm


[x,y,slice,freq, b1] = size(images);
freq_use = length(Para.Freq_ppm_use);
B0map = zeros(x,y,slice,1,b1)-3;
delta_freq_map = zeros(x,y,slice,1,b1)-3;
z_B0correction = zeros(x,y,slice,freq_use, b1);

for B1_loop = 1:b1
    for Slice_loop = Para.SelectedSlice
        img_S0 = medfilt2(images(:,:,Slice_loop,1,B1_loop));
        img_S0 = squeeze(img_S0);
        for f = 2:freq
        img_CEST(:,:,f-1) = medfilt2(images(:,:,Slice_loop,f,B1_loop));
        end      
        mask = Mask(:,:,Slice_loop);
            
            for ii = 1:x
                for jj = 1:y
                    if mask(ii, jj)==0
                        continue
                    end
                    B0_CEST = squeeze(img_CEST(ii,jj,:));
                    B0_S0 = squeeze(img_S0(ii,jj));
                    %   interpolate 2 -- 20230925
                    interp_zspectrum = interp1(Para.Freq_ppm,B0_CEST,Para.Freq_ppm_use,'cubic');  
                    
                    [B0map(ii,jj,Slice_loop,1,B1_loop),delta_freq_map(ii,jj,Slice_loop,1,B1_loop),zspectrum2fit] = B0_correction(interp_zspectrum',B0_S0,Para.Freq_ppm_use);
                    z_B0correction(ii,jj,Slice_loop,:,B1_loop) = zspectrum2fit;
                end
            end     
    end
end

end
