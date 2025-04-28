function [map_Cr,map_PCr,map_glycoNOE,res_B0map]=Fitting_PixelWise(images,Mask_all,B0map,Para,result_path)
%  calculate the GlycoNOE signals from the z_fit data
% including maps, signals summation, signal_area

% input:
% images: x,y,b1,freq, slice
% Mask_all:x,y,slice
% B0map: x,y,slice,1,b1
% B1map: x,y,slice

% output:
% map_Cr/map_PCr/map_glycoNOE: x,y,slice,b1
% Signal_Cr/Signal_PCr/Signal_glycoNOE: slice,b1

% xuxi updated 2023.8.4@siat 
% according to the codes from chongxue.bie

[nx, ny, nB1, nFreq, nSlice] = size(images);
images_Corr_z = zeros(size(images));
res_B0map = zeros(size(B0map))-3;

% initialize the maps
map_Cr = zeros(nx, ny,nSlice,nB1);
map_PCr = zeros(nx, ny,nSlice,nB1);
map_glycoNOE = zeros(nx, ny,nSlice,nB1);
% initialize the Signals
Signal_Cr = zeros(nSlice,nB1);
Signal_PCr = zeros(nSlice,nB1);
Signal_glycoNOE = zeros(nSlice,nB1);

% initialize the fitting images
images_downfield          = zeros(size(images));
images_upfield               = zeros(size(images)) ;
images_fit_back_down   = zeros(size(images)) ;
images_fit_back_up        = zeros(size(images)) ;

images_fit                        = zeros(size(images));
images_fit_p20                 = zeros(size(images));
images_fit_p25                 = zeros(size(images));
images_fit_no_p20           = zeros(size(images));
images_fit_no_p25           = zeros(size(images));
images_fit_no_Cr_PCr      = zeros(size(images));
images_fit_no_glycoNOE = zeros(size(images));
images_fit_up_mt            = zeros(size(images));
images_fit_down_mt       = zeros(size(images));

freq_ppm = Para.Freq_ppm;
sindex = [1:1:length(freq_ppm)];
n = 0;
%% calculate the CEST signal in each pixel
for B1 = 1:nB1
    for Slice = Para.SelectedSlice
        sum_Cr=0;% all the signals in ROI
        sum_PCr=0;
        sum_glycoNOE=0;
        
        mask = squeeze(Mask_all(:,:,Slice));
        
        for ii = 1:nx
            for jj = 1:ny
                if mask(ii, jj)==1                    
                    zspectrum2fit = squeeze(images(ii,jj,B1,:,Slice));
                    hasNaN = any(isnan(zspectrum2fit));
                    hasNegative = any(zspectrum2fit < 0);
                    if hasNaN | hasNegative
                         mask(ii, jj) = 0; 
                        continue
                    end
%                     zspectrum2fit(zspectrum2fit < 0 | isnan(zspectrum2fit)) = 0;

                    images_Corr_z(ii,jj,B1,:,Slice) = zspectrum2fit;% the original z spectrum after correction
                    res_B0map(ii,jj,Slice,1,B1)  = B0map(ii,jj,Slice,1,B1);
                    %% fit down-field Z-spectra, 5-pool: water, +2, +2.5, +3.5, MT
                    options = optimset('MaxFunEvals',100000*5, 'MaxIter',1000000,'Display','off');
                    
                    %           x = [  water_direct  |   +2                            |+2.5                         |+3.5                        |MT           dc ]
                    x  = [0.7       0        1     0.01     1.9      0.1       0.01     2.5    0.1       0.02     3.6        0.2      0.2     3.2    5      0.03];
                    lb = [0.05  -0.15     0       0        1.7      0         0          2.3    0           0        3.0          0        0       3.2    3     -0.05];
                    ub = [1      0.15      3     0.03     2.1      0.3       0.03     2.7    0.3       0.05     4.0        0.8      1.0      5     20     0.3];
                    pool_num = 5;
                    
                    fit_range_down_back = freq_ppm>0.2 & freq_ppm<=max(freq_ppm); % 0ppm is discarded
                    global signal_cal peak_cal
                    [x,fmin,residual] = lsqnonlin(@(x)multi_lorentz_fit_step1(x,freq_ppm(fit_range_down_back),zspectrum2fit(fit_range_down_back),pool_num),x,lb,ub,options);
                    
                    images_fit_back_down(ii,jj,B1,sindex(fit_range_down_back),Slice)  = signal_cal+peak_cal(:,2)+peak_cal(:,3)+peak_cal(:,4)+peak_cal(:,5)+x(end);%+peak_cal(:,6);
                    background_down = squeeze(images_fit_back_down(ii,jj,B1,:,Slice));
                    images_fit_down_mt(ii,jj,B1,sindex(fit_range_down_back),Slice) = peak_cal(:,5)+x(end);
                    
                    %% get down-field residual spectrum, and fit 4-pool:
                    downfield_spectrum = zeros(size(freq_ppm));                  
                    downfield_spectrum(freq_ppm>=0)=(background_down(freq_ppm>=0)-zspectrum2fit(freq_ppm>=0));
                    
                    images_downfield(ii,jj,B1,:,Slice) = downfield_spectrum; % sum of all the signals
                    
                    peak_cal = 0; signal_cal = 0;
                    %x = [          +2                              +2.5                        +3.5                        mT                 dc ]
                    x  =  [ 0.01    1.9      0.1       0.01     2.5    0.1       0.02     3.6        0.2    0.2     3.2    5        0];
                    lb =  [0         1.7      0            0        2.3    0         0          3.2        0      0        3.2    3      -0.05];
                    ub = [0.03     2.1      0.3       0.03     2.7    0.3       0.05     3.8        1.5    1.0     5    20      0.3];
                    pool_num = 4;
                    
                    fit_range_down = freq_ppm<=max(freq_ppm) & freq_ppm>0.3;
                    [x,fmin,residual] = lsqnonlin(@(x)multi_lorentz_fit_step2(x,freq_ppm(fit_range_down),downfield_spectrum(fit_range_down),pool_num),x,lb,ub,options);
                    
                    images_fit_p25(ii,jj,B1,sindex(fit_range_down),Slice) = peak_cal(:,2);
                    images_fit_p20(ii,jj,B1,sindex(fit_range_down),Slice) = peak_cal(:,1);
                    images_fit_no_p20(ii,jj,B1,sindex(fit_range_down),Slice) = peak_cal(:,2)+peak_cal(:,3)+x(end)+peak_cal(:,4);
                    images_fit_no_p25(ii,jj,B1,sindex(fit_range_down),Slice) = peak_cal(:,1)+peak_cal(:,3)+x(end)+peak_cal(:,4);
                    images_fit_no_Cr_PCr(ii,jj,B1,sindex(fit_range_down),Slice) = peak_cal(:,3)+x(end)+peak_cal(:,4);
                    sumindex_Cr = (freq_ppm<=2.5 & freq_ppm>=1.5);
                    Cr_signal = 100*(-squeeze(images_fit_no_p20(ii,jj,B1,sindex(sumindex_Cr),Slice))' + squeeze(downfield_spectrum(sumindex_Cr)));
                    
                    map_Cr(ii,jj,Slice,B1)=sum(Cr_signal);% all the signals in 1.5~2.5ppm
%                     sum_Cr = sum_Cr + sum(Cr_signal); % all the signals in ROI
                    
                    sumindex_PCr = (freq_ppm<=3 & freq_ppm>=2);
                    PCr_signal = 100*(-squeeze(images_fit_no_p25(ii,jj,B1,sumindex_PCr,Slice))' + squeeze(downfield_spectrum(sumindex_PCr)));
                    map_PCr(ii,jj,Slice,B1)=sum(PCr_signal);
%                     sum_PCr = sum_PCr + sum(PCr_signal);
                    
                    %% fit up-field Z-spectra, 4-pool: water, -1, -3, MT
                    peak_cal = 0; signal_cal = 0;
                    
                    %         x = [         water_direct,             -1                            -2.5                        -3.5                         MT           dc ]
                    x  = [0.7     0         1     0.02     -1        0.2       0.03     -2.5      0.2         0.02    -3.5     0.2         0.2     -3.2      5     0.03];
                    lb = [0.05   -0.15      0     0        -1.1      0         0        -2.9      0           0       -4.0     0           0       -5      3     -0.05];
                    ub = [1       0.15      3     0.05     -0.9      0.55      0.1      -2.0      0.6         0.2     -3.1     0.6         1.0     -3.2      20     0.3];
                    pool_num = 5;
                    
                    fit_range_up_back = freq_ppm>=min(freq_ppm) & freq_ppm<-0.2;
                    
                    [x,fmin,residual] = lsqnonlin(@(x)multi_lorentz_fit_step1(x,freq_ppm(fit_range_up_back),zspectrum2fit(fit_range_up_back),pool_num),x,lb,ub,options);
                    
                    images_fit_back_up(ii,jj,B1,sindex(fit_range_up_back),Slice) = signal_cal+peak_cal(:,2)+peak_cal(:,3)+peak_cal(:,4)+peak_cal(:,5)+x(end);%+peak_cal(:,6);                    
                    background_up =  squeeze(images_fit_back_up(ii,jj,B1,:,Slice));
                    images_fit_up_mt(ii,jj,B1,sindex(fit_range_up_back),Slice) = peak_cal(:,5)+x(end);
                    % xuxi updated
                    upfield_spectrum = zeros(size(freq_ppm));
                    upfield_spectrum(freq_ppm>0) = 0;
                    upfield_spectrum(freq_ppm<=0)=(background_up(freq_ppm<=0)-zspectrum2fit(freq_ppm<=0));                   
                    
                    images_upfield(ii,jj,B1,:,Slice) = upfield_spectrum;
                    
                    %%    up-field residual spectrum, and fit 3-pool:
                    
                    options = optimset('MaxFunEvals',100000*5, 'MaxIter',1000000,'Display','off');
                    % x = [          -1                                 -2.8                               -3.8                      MT           dc ]
                    x  = [  0.02     -1         0.2       0.03     -2.5      0.2         0.02    -3.5     0.2    0.2     -3.2      5     0];
                    lb = [  0          -1.1      0         0          -2.9       0           0         -4.0     0      0           -5     3    -0.05];
                    ub = [  0.05    -0.9      0.55      0.1      -2.0      0.6         0.2       -3.1     0.6    1.0     -3.2    20    0.3];
                    pool_num = 4;
                    
                    peak_cal = 0; signal_cal = 0;
                    fit_range_up = freq_ppm<-0.3 & freq_ppm>=min(freq_ppm);
                    [x,fmin,residual] = lsqnonlin(@(x)multi_lorentz_fit_step2(x,freq_ppm(fit_range_up),upfield_spectrum(fit_range_up),pool_num),x,lb,ub,options);
                   
                    images_fit_no_glycoNOE(ii,jj,B1,sindex(fit_range_up),Slice) =peak_cal(:,2)+peak_cal(:,3)+x(end) +peak_cal(:,4);
                    
                    sumindex_glycoNOE = (freq_ppm<=-0.6 & freq_ppm>=-1.6);
                    glycoNOE_signal = 100*(squeeze(images_upfield(ii,jj,B1,sindex(sumindex_glycoNOE),Slice)) - squeeze(images_fit_no_glycoNOE(ii,jj,B1,sindex(sumindex_glycoNOE),Slice)));
                    map_glycoNOE(ii,jj,Slice,B1)=sum(glycoNOE_signal);
%                     sum_glycoNOE = sum_glycoNOE + sum(glycoNOE_signal);
                    n = n +1;
                end
            end
        end
        
        
        %% ROI analysis
        
        cc = bwconncomp(mask,4); % find tubes
        tubes = regionprops(cc,'Centroid','Area'); % tube centers
        
        for nm = 1:numel(cc.PixelIdxList)
            tubes(nm).mask = zeros(nx, ny);
            tuberegion = cc.PixelIdxList{nm};
            tubes(nm).mask(tuberegion) = 1;
            Tarea = tubes(nm).Area;
            tube_mask = tubes(nm).mask;
            tubezspectrum(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_Corr_z(:,:,B1,:,Slice)))./Tarea);% the original z spectrum after correction
            tubezspectrum_resd_down(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_downfield(:,:,B1,:,Slice)))./Tarea); % sum of all the pools' signal 
            tubezspectrum_resd_up(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_upfield(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_back_down(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_back_down(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_back_up(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_back_up(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_p20(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_p20(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_p25(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_p25(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_no_p20(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_no_p20(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_no_p25(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_no_p25(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_no_Cr_PCr(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_no_Cr_PCr(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_no_glycoNOE(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_no_glycoNOE(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_up_mt(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_up_mt(:,:,B1,:,Slice)))./Tarea);
            tubezspectrum_fit_down_mt(nm,B1,:,Slice) = squeeze(sum(sum(tube_mask.*images_fit_down_mt(:,:,B1,:,Slice)))./Tarea);
            
            %% plot Z-spectrum
            
            fig = figure(1000+B1*100+Slice);
            subplot(4,2,[3,4,5,6]);
            hold on; box on;
            
            yyaxis left
            axis tight;
            hold on;
            plot(freq_ppm(2:end),squeeze(tubezspectrum(nm,B1,2:end,Slice)),'ko-','markersize',8,'LineWidth',1.5);
            hold  on
            
            plot(freq_ppm(fit_range_down_back),smooth(squeeze(tubezspectrum_fit_back_down(nm,B1,fit_range_down_back,Slice)),1),'b--', 'linewidth',2);
            plot(freq_ppm(fit_range_up_back),smooth(squeeze(tubezspectrum_fit_back_up(nm,B1,fit_range_up_back,Slice)),1),'b--', 'linewidth',2);
            
            set(gca,'fontsize',26)
            ylabel('S/S_0','fontsize',34);
            set(gca,'ylim', [0.4 1.05],'ytick',[0:0.2:1]);
            set(gca,'ycolor','k');
            set(gca,'ycolor','default')
            
            yyaxis right
            axis tight;
            hold on;
            box on;
            
            ylabel('Signal (%)','fontsize',34)
            xlabel('Offset from water (ppm)','fontsize',34)
            set(gca,'ycolor','b');
            
            plot(freq_ppm(fit_range_up_back),100*(squeeze(tubezspectrum_resd_up(nm,B1,(fit_range_up_back),Slice))),'bo-','markersize',8,'LineWidth',1.5);
            plot(freq_ppm(fit_range_up),100*(squeeze(tubezspectrum_fit_no_glycoNOE(nm,B1,(fit_range_up),Slice))),'r-', 'linewidth',2); 
            hold on
            
            a = freq_ppm(sumindex_glycoNOE);
            l1 = squeeze(tubezspectrum_resd_up(nm,B1,(sumindex_glycoNOE),Slice));
            l2 = squeeze(tubezspectrum_fit_no_glycoNOE(nm,B1,(sumindex_glycoNOE),Slice));
            fx = [a, fliplr(a)];
            fy = [100*l1', fliplr(100*l2')];
            h3 = fill(fx,fy,'r','marker','none');
            alpha(h3,0.3);
            
            plot(freq_ppm(fit_range_down_back),100*(squeeze(tubezspectrum_resd_down(nm,B1,(fit_range_down_back),Slice))),'bo-','markersize',8,'LineWidth',1.5);  
            plot(freq_ppm(fit_range_down),100*(squeeze(tubezspectrum_fit_no_Cr_PCr(nm,B1,(fit_range_down),Slice))),'r-', 'linewidth',2);  
            hold on
            
            a = freq_ppm(sumindex_Cr);
            l1 = squeeze(tubezspectrum_resd_down(nm,B1,(sumindex_Cr),Slice)); 
            l2 = squeeze(tubezspectrum_fit_no_Cr_PCr(nm,B1,(sumindex_Cr),Slice)); 
            fx = [a, fliplr(a)];
            fy = [100*l1', fliplr(100*l2')];
            h3 = fill(fx,fy,[0, 0.5, 0],'marker','none');
            alpha(h3,0.3);
            
            a = freq_ppm(sumindex_PCr); 
            l1 = squeeze(tubezspectrum_resd_down(nm,B1,(sumindex_PCr),Slice));
            l2 = squeeze(tubezspectrum_fit_no_Cr_PCr(nm,B1,(sumindex_PCr),Slice));
            fx = [a, fliplr(a)];
            fy = [100*l1', fliplr(100*l2')];
            h3 = fill(fx,fy,[0.75, 0, 0.75],'marker','none');
            alpha(h3,0.3);
            
            set(gca,'ylim', [0 20],'ytick',[0:5:80],'xlim',[-4 4],'xtick',[-8:1:8]);
            set(gca,'XDir','reverse')
            pos = get(gca, 'Position');
            pos(1) = 0.2;
            pos(3) = 0.6;
            set(gca, 'Position', pos)
            h=gcf;
            
            set(h,'Position',[800 300 900 800]);
            set(h,'PaperOrientation','landscape');
            saveas(h, [result_path,'z_B1_', num2str(B1),'_S_',num2str(Slice),'.png'], 'png');
            
            
            %%  spectrum imshow --important
            fig = figure( 2000+B1*100+Slice);
            subplot(4,2,[3,4,5,6]);
            hold on; box on;
            
%             zspectra_plot = squeeze(tubezspectrum(nm,1,sindex(indexf),Slice));
            yyaxis left
            axis tight;
            hold on;
            plot(freq_ppm(2:end),squeeze(tubezspectrum(nm,B1,2:end,Slice)),'ko-','markersize',8,'LineWidth',1.5);
            hold  on
            
            plot(freq_ppm(fit_range_down_back),smooth(squeeze(tubezspectrum_fit_back_down(nm,B1,fit_range_down_back,Slice)),1),'b--', 'linewidth',2);
            plot(freq_ppm(fit_range_up_back),smooth(squeeze(tubezspectrum_fit_back_up(nm,B1,fit_range_up_back,Slice)),1),'b--', 'linewidth',2);
            
            set(gca,'fontsize',26)
            ylabel('S/S_0','fontsize',34);
            set(gca,'ylim', [0.3 1.05],'ytick',[0:0.2:1]);
            set(gca,'ycolor','k');
            set(gca,'ycolor','default')
            
            yyaxis right
            axis tight;
            hold on;
            box on;
            
            ylabel('Signal (%)','fontsize',34)
            %        xlabel('Offset from water (ppm)','fontsize',34)
            set(gca,'ycolor','b');
            %        set(gca,'ylim', [-2 8],'ytick',[0:5:15]);
%             baseline = zeros(length(freq_ppm));
            %        plot(freq_ppm,baseline,'k--', 'linewidth',2);
            
            plot(freq_ppm(fit_range_up_back),100*(squeeze(tubezspectrum_resd_up(nm,B1,(fit_range_up_back),Slice))),'marker','o','markersize',8,'linestyle','-','color','b', 'linewidth',1.5);
            plot(freq_ppm(fit_range_up),100*(squeeze(tubezspectrum_fit_no_glycoNOE(nm,B1,(fit_range_up),Slice))),'marker','none','linestyle','-','color','r', 'linewidth',2);
            hold on
            
            plot(freq_ppm(fit_range_down_back),100*(squeeze(tubezspectrum_resd_down(nm,B1,fit_range_down_back,Slice))),'marker','o','markersize',8,'linestyle','-','color','b', 'linewidth',1.5);
            plot(freq_ppm(fit_range_down),100*(squeeze(tubezspectrum_fit_no_Cr_PCr(nm,B1,fit_range_down,Slice))),'marker','none','linestyle','-','color','r', 'linewidth',2);
            
            
            set(gca,'ylim', [0 20],'ytick',[0:5:80],'xlim',[-4 4],'xtick',[-8:1:8]);
            set(gca,'XDir','reverse')
            %        set(gca,'fontsize',21)
            pos = get(gca, 'Position');
            pos(1) = 0.2;
            pos(3) = 0.6;
            set(gca, 'Position', pos)
            %        set(gca,'fontsize',25)
            
            subplot(subplot(4,2,7:8)); hold on; box on;
            pos = get(gca, 'Position');
            pos(1) = 0.2;
            pos(3) = 0.6;
            set(gca, 'Position', pos)
            
            %      glycoNOE
            plot(freq_ppm(fit_range_up)',100*(smooth(squeeze(tubezspectrum_resd_up(nm,B1,(fit_range_up),Slice)),1)-smooth(squeeze(tubezspectrum_fit_no_glycoNOE(nm,B1,(fit_range_up),Slice)),1)),'r-', 'linewidth',2);
            sum_glycoNOE = 100*(smooth(squeeze(tubezspectrum_resd_up(nm,B1,(sumindex_glycoNOE),Slice)),1)-smooth(squeeze(tubezspectrum_fit_no_glycoNOE(nm,B1,(sumindex_glycoNOE),Slice)),1));
            glycoNOE_area =sum (sum_glycoNOE) % 前面算了一个是干嘛的？
            Signal_glycoNOE(Slice,B1) = glycoNOE_area; % area in the ppm range in ROI
            
            area_plot = [ sum_glycoNOE];
            h(Slice) = area(freq_ppm(sumindex_glycoNOE), area_plot);
            h(Slice).FaceColor = 'r';
            h(Slice).EdgeColor = 'r';
            alpha(h(Slice),0.3)
            hold on
            
            peak_glycoNOE = (freq_ppm<=-0.99 & freq_ppm>=-1.11);
            glycogen_spec = 100*(smooth(squeeze(tubezspectrum_resd_up(nm,B1,:,Slice)),1) - smooth(squeeze(tubezspectrum_fit_no_glycoNOE(nm,B1,:,Slice)),1));
            glycogen_amplitude(Slice,B1) = mean(glycogen_spec(peak_glycoNOE));
            
            % Cr
            fit_range_Cr = freq_ppm<4.9 & freq_ppm>0.9;
            plot(freq_ppm(fit_range_Cr),100*(-smooth(squeeze(tubezspectrum_fit_no_p20(nm,B1,fit_range_Cr,Slice)),1) + smooth(squeeze(tubezspectrum_resd_down(nm,B1,fit_range_Cr,Slice)),1)),'linestyle','-','color',[0 0.5 0], 'linewidth',2);
            
            sum_Cr = 100*(-smooth(squeeze(tubezspectrum_fit_no_p20(nm,B1,sumindex_Cr,Slice)),1) + smooth(squeeze(tubezspectrum_resd_down(nm,B1,sumindex_Cr,Slice)),1));
            Cr_area =sum (sum_Cr) 
            Signal_Cr(Slice,B1) = Cr_area;
            
            area_plot = [ sum_Cr];
            h(Slice+1) = area(freq_ppm(sumindex_Cr), area_plot);
            h(Slice+1).FaceColor = [0, 0.5, 0];
            h(Slice+1).EdgeColor = [0, 0.5, 0];
            alpha(h(Slice+1),0.3)
            peak_Cr = (freq_ppm<=2.01 & freq_ppm>=1.89);
            Cr_spec = 100*(-smooth(squeeze(tubezspectrum_fit_no_p20(nm,B1,:,Slice)),1) + smooth(squeeze(tubezspectrum_resd_down(nm,B1,:,Slice)),1));
            Cr_amplitude(Slice,B1) = mean(Cr_spec(peak_Cr)); 
            
            
            %      PCr
            fit_range_PCr = freq_ppm<4.9 & freq_ppm>0.9;
            plot(freq_ppm(fit_range_PCr),100*(-smooth(squeeze(tubezspectrum_fit_no_p25(nm,B1,fit_range_PCr,Slice)),1) + smooth(squeeze(tubezspectrum_resd_down(nm,B1,fit_range_PCr,Slice)),1)),'marker','none','linestyle','-','color',[0.75 0 0.75], 'linewidth',2);
            sum_PCr = 100*(-smooth(squeeze(tubezspectrum_fit_no_p25(nm,B1,sumindex_PCr,Slice)),1) + smooth(squeeze(tubezspectrum_resd_down(nm,B1,sumindex_PCr,Slice)),1));
            PCr_area =sum (sum_PCr) 
            Signal_PCr(Slice,B1) = PCr_area;
            
            area_plot = [ sum_PCr];
            h(Slice+2) = area(freq_ppm(sumindex_PCr), area_plot);
            h(Slice+2).FaceColor = [0.75, 0, 0.75];
            h(Slice+2).EdgeColor = [0.75, 0, 0.75];
            alpha(h(Slice+2),0.3)
            
            peak_PCr = (freq_ppm<=2.61 & freq_ppm>=2.49);
            PCr_spec = 100*(-smooth(squeeze(tubezspectrum_fit_no_p25(nm,B1,:,Slice)),1) + smooth(squeeze(tubezspectrum_resd_down(nm,B1,:,Slice)),1));
            PCr_amplitude(Slice,B1) = mean(PCr_spec(peak_PCr)); 
            
            
            set(gca,'fontsize',26)
            set(gca,'XDir','reverse')
            set(gca,'ylim', [-0.5 4],'ytick',[-5:1:25],'xlim',[-4 4],'xtick',[-6:1:6]);
            xlabel('Offset from water (ppm)','fontsize',34)
            ylabel('Signal (%)','fontsize',34)
            grid on;
            h=gcf;
            
            set(h,'Position',[800 300 1000 1300]);
            set(h,'PaperOrientation','landscape');
            saveas(h, [result_path,'signal_B1_', num2str(B1),'_S_',num2str(Slice),'.png'], 'png');
        end
    end
end

%% save the signals and amplitdes of CEST
save([result_path,'Signal.mat'],'Signal_Cr','Signal_PCr','Signal_glycoNOE'); % [Slice,B1]

save([result_path,'Map.mat'],'map_Cr','map_PCr','map_glycoNOE','res_B0map'); % [x,y,Slice,B1]

save([result_path,'Amp.mat'],'glycogen_amplitude','Cr_amplitude', 'PCr_amplitude','tubezspectrum'); % [Slice,B1],  [B1, freq,Slice] = size(tubezspectrum);
end