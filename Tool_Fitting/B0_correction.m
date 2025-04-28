function [B0map,delta_freq_map,zspectrum2fit] = B0_correction(cestData,S0Data,freq_ppm)
% cestData:[freq]
% S0Data:[1]
% delta_freq_map: 0ppm offset
% B0map: delta_B0 value
% zspectrum2fit: z spectrum after B0 correction
    sindex = [1:1:length(freq_ppm)];
    int = squeeze(1-cestData(sindex)/max(cestData(sindex)));
    para(2)=0;
    
    gyr = 42.58 *10^6;
    fit_freq_range = (freq_ppm > -1.5 & freq_ppm < 1.5); 
    interp_ratio = 5;  
    
    interp_int = interp(int(fit_freq_range), interp_ratio); 
    interp_fit_freq_ppm = interp(freq_ppm(fit_freq_range),interp_ratio);
    [fit,para] = lorentzfit(interp_fit_freq_ppm,interp_int',[0, 0, 1, 0.1]);
    delta_freq = para(2); % get zero-value frequency
    delta_freq_map = delta_freq; 
    B0map = -para(2)/gyr; % tranform to B0 value
%     zspectrum = squeeze(cestData(sindex(1:end)))';
            
    interp_zspectrum = interp(cestData',interp_ratio);
    interp_freq_ppm = interp(freq_ppm,interp_ratio);
    interp_zspectrum_corrected = zspectrum_corr(interp_freq_ppm,-delta_freq,interp_zspectrum);
    zspectrum_corrected = downsample(interp_zspectrum_corrected ,interp_ratio)'; 
            
    zspectrum2fit = zspectrum_corrected/S0Data;
end


function zspectrum_corrected = zspectrum_corr(freq,delta_freq,zspectrum)
zspectrum_corrected = zeros(length(freq),1);
for ii = 1:length(freq)    
    [~,index] = min((freq + delta_freq  - freq(ii)).^2);% the closest point after translation by 0 ppm offset
    zspectrum_corrected(ii) = zspectrum(index);
end

end