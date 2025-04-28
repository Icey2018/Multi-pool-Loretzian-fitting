function ShowSave_map(map_Cr,map_PCr,map_glycoNOE,B0map,B1map,T2w,images, Mask,result_path,Para,images_registraiton,Mask1)
% show the maps in ROI, if there isn't T2w, show S0; ones(x,y) in case of no  B1 map

% map_Cr: x,y,slice,b1
% res_B0map:x,y,slice,1,b1
% B1map_registraiton: x,y,slice
% images: x,y,slice, freq, b1
% Mask: x,y,slice


% T2w_registraiton: [x,y,slice]


[x,y,Slice,B1] = size(map_Cr);

for b =1:B1
    for s = Para.SelectedSlice
        fig = figure( 3000+b*100+s);
        % T2
        ax1 = subplot(4,2,1);
        
        if ~isempty(T2w)
            imagesc(imresize(T2w(:,:,s),[x,y]));
            %imshow(imresize(T2w(:,:,s), [x, y]));%wuhao
        else
            imagesc(images(:,:,b,1,s)); %  images: x,y,slice, freq, b1
            %imshow(images(:,:,b,1,s));%wuhao
        end
        %        colormap(ax1,gray)
        altered_colormap(ax1, gray) % test %进行色彩映射调整
        ax1.Position = ax1.Position + [0 0 -0.065 0];
        
        title('T_2','FontSize',22);
        axis off;
        hold on
        contour(gca, squeeze(Mask(:,:,s)), [0.5 0.5],'r-','LineWidth',2);%画线

        %{
        imshow(squeeze(images_registraiton(:,:,s,1,b)),[]);
        altered_colormap(ax1, gray)
        ax1.Position = ax1.Position + [0 0 -0.065 0];
        title('T_2','FontSize',22);
        axis off;
        hold on
        contour(gca, squeeze(Mask(:,:,s)), [0.5 0.5],'r-','LineWidth',2);

        %}
        % B1 map
        ax2 = subplot(4,2,2);
        if ~isempty(B1map)
            imagesc(B1map(:,:,s)*100.*squeeze(Mask(:,:,s)));
        else
            imagesc(ones(x,y));
        end
        altered_colormap(ax2,jet(256));
        caxis([60 140]);
        ax2.Position = ax2.Position + [0 0 0.01 0];
        axis off
        pos = get(gca, 'Position');
        set(gca, 'Position', pos)
        h = colorbar( ...
            'FontSize',20);
        title('B_1 map (%)','FontSize',22);
        
        % B0 map
        ax3 = subplot(4,2,3);
        B0map_show=medfilt2(squeeze(B0map(:,:,s,1,b)),[4,4]);
        imagesc(B0map_show);
        colormap(ax3, jet(256))
        altered_colormap(ax3,jet(256))
        caxis([-1 1]);
        axis off
        pos = get(gca, 'Position');
        set(gca, 'Position', pos)
        h = colorbar( ...
            'FontSize',20);
        title('B_0 map (ppm)','FontSize',22);
        
        % Cr map
        ax4 = subplot(4,2,4);
        map_Cr_show=medfilt2(squeeze(map_Cr(:,:,s,b)),[4,4]);
        imagesc(map_Cr_show);
        altered_colormap(ax4,jet(256));
        axis off
        caxis([0 25]);
        pos = get(gca, 'Position');
        set(gca, 'Position', pos)
        h = colorbar( ...
            'FontSize',20);
        set(get(h,'ylabel'),'string','Contrast (%*ppm)','FontSize',20);
        title('+1.95 ppm','FontSize',22);
        % PCr map
        ax5 = subplot(4,2,5);
        map_PCr_show=medfilt2(squeeze(map_PCr(:,:,s,b)),[4,4]); 
        imagesc(map_PCr_show);
        altered_colormap(ax5,jet(256))
        axis off
        caxis([0 25]);
        pos = get(gca, 'Position');
        set(gca, 'Position', pos)
        h = colorbar( ...
            'FontSize',20);
        set(get(h,'ylabel'),'string','Contrast (%*ppm)','FontSize',20);
        title('+2.5 ppm','FontSize',22);
        
        % glycoNOE map
        ax6 = subplot(4,2,6);
        map_glycoNOE_show=medfilt2(squeeze(map_glycoNOE(:,:,s,b)),[4,4]);
        imagesc(map_glycoNOE_show);
        altered_colormap(ax6,jet(256))
        axis off
        caxis([0 40]);%设置当前子图的色彩轴范围为从 0 到 40
        pos = get(gca, 'Position');
        set(gca, 'Position', pos) %获取当前子图的坐标轴位置，并将其设置回相同的位置
        h = colorbar( ...
            'FontSize',20);
        set(get(h,'ylabel'),'string','Contrast (%*ppm)','FontSize',22);
        title('glycoNOE','FontSize',22);

        ax1 = subplot(4,2,7);
        imshow(squeeze(images_registraiton(:,:,s,1,b)),[]);
        altered_colormap(ax1, gray)
        ax1.Position = ax1.Position + [0 0 -0.065 0];
        axis off;
        hold on
        contour(gca, squeeze(Mask(:,:,s)), [0.5 0.5],'r-','LineWidth',2);
        %contour(gca, squeeze(Mask1(:,:,s)), [0.5 0.5],'g-','LineWidth',2);
        pos = get(gca, 'Position');
        set(gca, 'Position', pos)
        %h = colorbar( ...
            %'FontSize',20);
        title('ROI','FontSize',22);
        
        % figure
        h=gcf;
        set(h,'Position',[500 50 600 1000]);
        set(h,'PaperOrientation','landscape');
        saveas(h, [result_path,'Map_B1_',num2str(b),'_S_',num2str(s),'.png'], 'png');
    end
end

end