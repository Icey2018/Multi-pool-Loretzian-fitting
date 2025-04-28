function [f] = multi_lorentz_fit_step1(x, freq_ppm, zspectrum2fit, n)
% output:
% signal_cal : 1-减去所有池和dc的剩余信号
% peak_cal：各个池的信号幅值
% f: z_fit-z spectrum, the residual 

% input:
%  x: para of fitting
%  zspectrum2fit: z spectrum value
%  n: pool number
%  interp_ratio: default value is 1,  if the points number is not enough, the interp_ratio need to be bigger


global signal_cal peak_cal


a  = x(1:3:3*n); % amplitude Ai
w0 = x(2:3:3*n);% chemical shift wi
h  = x(3:3:3*n); % linewidth sigma
dc = x(end);% background (?)

signal_cal = ones(size(zspectrum2fit))-dc; % 1-dc
for ii = 1:numel(freq_ppm)
    for jj = 1:n
        peak_cal(ii,jj) = a(jj)/( 1 + (freq_ppm(ii) - w0(jj))^2/(0.5* h(jj))^2 ); % equation
        signal_cal(ii) = signal_cal(ii) - peak_cal(ii,jj);   % 1-dc-all_pool
    end
end

f = signal_cal - zspectrum2fit; % 1-dc-sum(pool) = z-spectrum（definition）

end

