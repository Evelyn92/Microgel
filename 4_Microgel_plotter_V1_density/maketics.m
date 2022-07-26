function [anything]=maketics(ticdist)
% Eric
%set(gca,'Xtick',[0:ticdist:1000]);
%set(gca,'XtickLabel',[0:ticdist:1000]);
% set(gca,'Ytick',[0:ticdist:1000]);
% set(gca,'YtickLabel',[0:ticdist:1000]);
% set(gca,'Ztick',[0:ticdist:1000]);
% set(gca,'ZtickLabel',[0:ticdist:1000]);

% Ashvini Modified
% set(gca,'XtickLabel',{-300:ticdist:300});
set(gca,'XtickLabel',{'-300';'-150';'0';'150';'300'});%To put it right below the tick
set(gca,'Xtick',[-300:ticdist:300],'TickDir','out', 'TickLength', [0.02 0.035]);

set(gca,'Ytick',[-300:ticdist:300],'TickDir','out', 'TickLength', [0.02 0.035]);
%set(gca,'YtickLabel',[-300:ticdist:300]);
set(gca,'YtickLabel',{'-300';'-150';'0';'150';'300'});%To put it right below the tick

set(gca,'Ztick',[-300:ticdist:300],'TickDir','out', 'TickLength', [0.02 0.035]);
set(gca,'ZtickLabel',[-300:ticdist:300]);
end