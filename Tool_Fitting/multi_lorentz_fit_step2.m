function [ f ] = multi_lorentz_fit_step2(x, freq_ppm, curve2fit, n)
global signal_cal peak_cal
% 整体思路：拟合后的信号signal_cal和输入进来的信号curve2fit逼近
% output:
% f: z_fit-z spectrum, the residual
% signal_cal: summary of all the pools
% peak_cal: each pool

% input:
%  x: para of fitting
%  zspectrum2fit: z spectrum value
%  n: pool number


a  = x(1:3:3*n);
w0 = x(2:3:3*n);
h  = x(3:3:3*n);
dc = x(end);

signal_cal = zeros(size(curve2fit))+dc;
for ii = 1:numel(freq_ppm)
    for jj = 1:n
        peak_cal(ii,jj) = a(jj)/( 1 + (freq_ppm(ii) - w0(jj))^2/(0.5* h(jj))^2 ); 
        signal_cal(ii) = signal_cal(ii) + peak_cal(ii,jj);% different from step1, this is the sum of CEST signals
    end
end

f = signal_cal - curve2fit;

end

