function [images_corr_motion, images]   = Registration_img(Para,images)
% registrated each images in each slice to ones beyond +-0.6ppm of the 1st
% B1 value
[x,y,nSlice,nFreq,nb1] = size(images);
freq_ppm = Para.Freq_ppm;
mtype = 'ssd';%'cc';
ttype = 'rigid';
depth = 4;

for loop = Para.SelectedSlice
    %         ref: averaged image between -0.6 to -1.6 ppm
    ref_ind = find((freq_ppm<=-0.6 & freq_ppm>=-1.6)|(freq_ppm>=0.6 & freq_ppm<=1.6));
    ref_ind = ref_ind + 1; %1st is S0.
    Vo_fix = mean(squeeze(images(:,:,loop,ref_ind,1)),3);
    %         do not register images close 0 ppm
    no_reg_ind = find(freq_ppm<=0.6 & freq_ppm>=-0.6);
    no_reg_ind = no_reg_ind + 1;
    
    %         reg saturated image to -1 ppm £¨mean(-0.6~-1.6ppm) £©
    for b = 1:nb1
        for ii = 1:length(freq_ppm)
            if ((ii+1) ~= no_reg_ind)
                Vo_moving = squeeze(images(:,:,loop,ii+1,b));
                tmp = image_registration_pairwise(Vo_moving,Vo_fix,mtype,ttype,depth);
                images_corr_motion(:,:,loop,ii+1,b) = tmp;
            else
                images_corr_motion(:,:,loop,ii+1,b) = images(:,:,loop,ii+1,b);
            end
        end        
        
        % reg S0 image to -1 ppm
        Vo_moving = squeeze(images(:,:,loop,1,b));
        tmp = image_registration_pairwise(Vo_moving,Vo_fix,mtype,ttype,depth);
        images_corr_motion(:,:,loop,1,b) = tmp;
    end
end
end