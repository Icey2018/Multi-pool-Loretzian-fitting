
function Mask = DrawMask(image,ROI_number)
color_all = [112,48,160;255,192,0;...
    91,155,213;146,208,80;
    192,0,0;255,230,153;214,143,233;...
    63,67,229]/255;

titleMessage = 'Choose ROI';
figure;imshow(abs(image),[],'InitialMagnification','fit');% show M0 image
for mm = 1:ROI_number % loop over different ROIs
    title([strrep(titleMessage,'_',' '),' ',num2str(mm)]);drawnow
    eval(['[Mask_',num2str(mm),',xindex_',num2str(mm),',yindex_',num2str(mm),'] = roipoly;']);% mutually choose
    eval(['Mask(:,:,',num2str(mm),')=Mask_',num2str(mm),';'])
    eval(['line(xindex_',num2str(mm),',yindex_',num2str(mm),',''Color'',color_all(',num2str(mm),',:),''LineWidth'',3.5);drawnow'])
end

end