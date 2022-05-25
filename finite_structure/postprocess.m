function postprocess(filename,width,height,fontsize)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Postprocesses figure to have right size and fontsize %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure will be saved in the same folder as postprocess
% width and height need to be specified in inches
% fontsize as usual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load file
[pathstr, name, ext] = fileparts(filename);
open(filename)
% Change size
set(gcf,'units','centimeters')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])
set(gcf, 'renderer', 'painters');
% Change font
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
% Export file
print('-depsc2',name)

end
