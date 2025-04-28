function altered_colormap(ax,cmap)
% cmap = colormap;
nzeros = 2;
cmap(1:nzeros,:) = zeros(nzeros,3);
colormap(ax,cmap); 