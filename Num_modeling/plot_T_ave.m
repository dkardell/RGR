clear

load T_ave_syn1.mat; T_ave_syn1 = T_ave;
load T_ave_syn2.mat; T_ave_syn2 = T_ave;
load T_ave_syn3.mat; T_ave_syn3 = T_ave;
load T_ave_syn4.mat; T_ave_syn4 = T_ave;
load T_ave_syn5.mat; T_ave_syn5 = T_ave;
load T_ave_crest.mat; T_ave_crest = T_ave;
clear T_ave

figure
set(gcf,'units','normalized','outerposition',[0 0 0.3 0.3])
plot(T_ave_syn1,'r','linewidth',2); hold on
plot(T_ave_syn2,'b','linewidth',2)
plot(T_ave_syn3,'g','linewidth',2)
plot(T_ave_syn4,'y','linewidth',2)
plot(T_ave_syn5,'m','linewidth',2)
plot(T_ave_crest,'k','linewidth',2); hold off
ylabel('Temperature [^oC]')
ylim([42 49])
legend('Syn1','Syn2','Syn3','Syn4','Syn5','CREST','Location','southeast')
% set(gca,'xtick',[],'ytick',[])
% print('T_eq','-dpsc')