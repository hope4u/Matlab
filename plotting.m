%% Generate Plots for MIMO Sem Project Marc M�ller
close all
%% cdfs of the achievavle Rate for wf-Solutions
load('MarcMIMO_VBLAST.combined.mat')
plots = {'none','wf','sp_iwf_paper','numericalGrad_VBlast'};
names = {'none','waterfilling','iterative waterfilling by jindal','numerical Gradient for max sumRate'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1;1 1 0];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('optimizing sumRate of V-BLAST at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_sum(snr,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
end

hold off
grid on
switch SNR(snr)
    case -10
        xlim([8,16]);
    case 0
        xlim([28,40]);
    case 20
        xlim([95,115]);
    case 40
        xlim([170,195]);
end
legend(leg,'Location','Best');
xlabel('sumRate');
ylabel('empirical CDF');
name = sprintf('results_vblast_none_wf_iwf_ngrad_%idb',SNR(snr));
saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
saveas(gcf,sprintf('figures/%s.fig',name))
% print('-depsc' myplot.eps
end

%% compare vblast and mmse
plots = {'none','wf','numericalGrad_VBlast','none','numericalGrad'};
names = {'none for V_BLAST','waterfilling','max sumRate for V-BLAST','none for LMMSE','max sumRate for LMMSE'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1;1 1 0];
none = 0;
for snr=1:numel(SNR)
figure
hold on
title(sprintf('Compare V-BLAST and LMMSE at %i dB SNR',SNR(snr)));
leg = {};
load('MarcMIMO_VBLAST.combined.mat')
for num=1:3
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_sum(snr,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
end
load('MarcMIMO_LMMSE.combined.mat')
for num=4:5
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_sum(snr,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
end

hold off
grid on
switch SNR(snr)
    case -10
        xlim([7,16]);
    case 0
        xlim([15,40]);
    case 20
        xlim([35,115]);
    case 40
        xlim([40,195]);
end
legend(leg,'Location','Best');
xlabel('sumRate');
ylabel('empirical CDF');
name = sprintf('results_vblast_mmse_none_wf_ngrad_%idb_2',SNR(snr));
saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
saveas(gcf,sprintf('figures/%s.fig',name))
end
%% cdfs for throughput Maximization V-Blast
load('MarcMIMO_VBLAST.combined.mat')
plots = {'none','numericalGrad_VBlast','minmax_VBlast'};
names = {'none','max sumRate','max minRate'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('maxmin vs max sumRate for V-BLAST at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_sum(snr,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
end

hold off
grid on
switch SNR(snr)
    case -10
        xlim([7,16]);
    case 0
        xlim([27,40]);
    case 20
        xlim([84,115]);
    case 40
        xlim([145,195]);
end
legend(leg,'Location','Best');
xlabel('sumRate');
ylabel('empirical CDF');
name = sprintf('results_vblast_none_ngrad_maxmin_%idb',SNR(snr));
saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
saveas(gcf,sprintf('figures/%s.fig',name))
end

%% cdfs for Rates V-BLAST
load('MarcMIMO_VBLAST.combined.mat')
plots = {'numericalGrad_VBlast','minmax_VBlast'};
names = {'none','max sumRate','max minRate'};
color = [1 0 1; 1 0 .5; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 1 0];
for snr=1:numel(SNR)
figure
hold on
title(sprintf('minimum/maximum Rates at %i dB SNR',SNR(snr)));
leg = {};

% for num=1:numel(plots)
num = 1;
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        sor = sort(squeeze(R_ch(:,snr,Opt,:)),1);
        for ch=1:12
            Vector = sor(ch,:);
            [f,x] = ecdf(Vector);
            plot(x,f,'-','LineWidth',3,'Color',color(ch,:));
            set(gca, 'Yscale', 'log');
            %leg{end+1}=sprintf('per channel Rates of %s',names{num});
        end
        leg{end+1}=sprintf('minRate of %s',names{num});
    end
% end
num = 2;
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(min(R_ch(:,snr,Opt,:)));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',[.3,.3,1]);
        leg{end+1}=sprintf('minRate of %s',names{num});

        Vector = squeeze(max(R_ch(:,snr,Opt,:)));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',[.3,.3,1]);
        set(gca, 'Yscale', 'log');
        leg{end+1}=sprintf('maxRate of %s',names{num});
    end
hold off
grid on
switch SNR(snr)
    case -10
        xlim([0,2]);
    case 0
        xlim([0,4]);
    case 20
        xlim([0,11]);
    case 40
        xlim([4,18]);
end
% legend(leg,'Location','Best');
xlabel('Rate');
ylabel('empirical CDF');
name = sprintf('results_vblast_ngrad_perCh_%idb',SNR(snr));
saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
saveas(gcf,sprintf('figures/%s.fig',name))
end
%% cdfs for throughput Maximization
load('MarcMIMO_LMMSE.combined.mat')
plots = {'none','numericalGrad','minmax'};
names = {'none','max sumRate','max minRate'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('maxmin vs max sumRate for LMMSE at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_sum(snr,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
end

hold off
grid on
switch SNR(snr)
    case -10
        xlim([7,12]);
    case 0
        xlim([15,32]);
    case 20
        xlim([35,95]);
    case 40
        xlim([40,180]);
end
legend(leg,'Location','Best');
xlabel('sumRate');
ylabel('empirical CDF');
name = sprintf('results_mmse_none_ngrad_maxmin_%idb_2',SNR(snr));
saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
saveas(gcf,sprintf('figures/%s.fig',name))
end
%% cdfs for Rates
load('MarcMIMO_LMMSE.combined.mat')
plots = {'numericalGrad','minmax'};
names = {'none','max sumRate','max minRate'};
color = [1 0 1; 1 0 .5; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 1 0];
for snr=1:numel(SNR)
figure
hold on
title(sprintf('minimum/maximum Rates at %i dB SNR',SNR(snr)));
leg = {};

% for num=1:numel(plots)
num = 1;
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        sor = sort(squeeze(R_ch(:,snr,Opt,:)),1);
        for ch=1:12
            Vector = sor(ch,:);
            [f,x] = ecdf(Vector);
            plot(x,f,'-','LineWidth',3,'Color',color(ch,:));
            set(gca, 'Yscale', 'log');
            %leg{end+1}=sprintf('per channel Rates of %s',names{num});
        end
        leg{end+1}=sprintf('minRate of %s',names{num});
    end
% end
num = 2;
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(min(R_ch(:,snr,Opt,:)));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',[.3,.3,1]);
        leg{end+1}=sprintf('minRate of %s',names{num});

        Vector = squeeze(max(R_ch(:,snr,Opt,:)));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',[.3,.3,1]);
        set(gca, 'Yscale', 'log');
        leg{end+1}=sprintf('maxRate of %s',names{num});
    end
hold off
grid on
% legend(leg,'Location','Best');
xlabel('Rate');
ylabel('empirical CDF');
name = sprintf('results_mmse_ngrad_perCh_%idb_2',SNR(snr));
saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
saveas(gcf,sprintf('figures/%s.fig',name))
end

% %% cdfs for power Minimization
% load('MarcMIMO_LMMSE.combined.mat')
% plots = {'powermin','fodorPrecoding'};
% names = {'numerical Gradient','Algorithm by Fodor'};
% Rates = [1,2,10];
% color = [1 0 0; 0 1 0; 0 0 1];
% 
% for snr=1:numel(SNR)
% figure
% hold on
% 
% title(sprintf('total Power for Fixed Target SINR at %i noise Power',1/10^(SNR(snr)/10)));
% leg = {};
% 
% for num=1:numel(plots)
%     Opt = find(ismember(Optimizer,plots(num)));
%     if sum(Opt)
%         for pl = 1:length(Opt)
% %             Vector = squeeze(R_sum(snr,Typ,Opt(pl),:));
% %             [f,x] = ecdf(Vector);
% %             plot(x,f,'-','LineWidth',3,'Color',color(num,:));
% 
%             Vector = squeeze(P_tot(snr,Opt(pl),:));
%             [f,x] = ecdf(Vector);
%             plot(x,f,'-','LineWidth',3,'Color',color(num,:)/pl);
%             set(gca, 'Yscale', 'log');
%             leg{end+1}=sprintf('%s at a Rate of %i dB',names{num},Rates(pl));
%         end
%     end
% end
% 
% hold off
% grid on
% legend(leg,'Location','Best');
% xlabel('Power');
% ylabel('empirical CDF');
% name = sprintf('results_mmse_powermin_all_%idb',SNR(snr));
% saveas(gcf,sprintf('figures/%s.eps',name))
% saveas(gcf,sprintf('figures/%s.fig',name))
% end
%% cdfs for power Minimization
load('MarcMIMO_LMMSE.combined.mat')
plots = {'powermin','fodorPrecoding'};
names = {'numerical Gradient','Algorithm by Fodor'};
Rates = [1,2,10];
color = [1 0 0; 0 1 0; 0 0 1];



for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    for pl = 1:length(Opt)
        figure
        hold on

        title(sprintf('total Power for Target Rate of %i', Rates(pl)));
        leg = {};
        for snr=1:numel(SNR)
            if snr>3
                continue
            end
%             Vector = squeeze(R_sum(snr,Typ,Opt(pl),:));
%             [f,x] = ecdf(Vector);
%             plot(x,f,'-','LineWidth',3,'Color',color(num,:));

            Vector = squeeze(P_tot(snr,Opt(pl),:));
            [f,x] = ecdf(Vector);
            plot(x,f,'-','LineWidth',3,'Color',color(num,:)/snr);
            set(gca, 'Yscale', 'log');
            leg{end+1}=sprintf('noise power %3.1f',1/10^(SNR(snr)/10));
        end
        hold off
        grid on
        switch Rates(pl)
            case 1
                xlim([0,30]);
            case 2
                xlim([0,200]);
            case 10
                xlim([0,3E5]);
        end
        legend(leg,'Location','Best');
        xlabel('Power');
        ylabel('empirical CDF');
        name = sprintf('results_mmse_%s_%i',plots{num},Rates(pl));
        saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
        saveas(gcf,sprintf('figures/%s.fig',name))
    end
end
%% ratios for power Minimization
load('MarcMIMO_LMMSE.combined.mat')
plots = {'powermin','fodorPrecoding'};
names = {'numerical Gradient','Algorithm by Fodor'};
Rates = [1,2,10];
color = [1 0 0; 0 1 0; 0 0 1];

Opt = find(ismember(Optimizer,plots(1)));
for pl = 1:length(Opt)
    
    figure
    hold on
    title(sprintf('total power ratios for target rate of %i', Rates(pl)))
    leg = {};

    for snr=1:numel(SNR)
        if snr>2
            continue
        end
        Values = P_tot(snr,Opt(pl)+length(Opt),:)./P_tot(snr,Opt(pl),:);
        Vector = squeeze(Values);
        
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(snr,:));
%         set(gca, 'Xscale', 'log');
        leg{end+1}=sprintf('noise power %3.1f',1/10^(SNR(snr)/10));
        
    end
    hold off
    grid on
        switch Rates(pl)
            case 1
                xlim([.5,1.5]);
            case 2
                xlim([0.9,1.1]);
            case 10
                xlim([0,100]);
        end
    legend(leg,'Location','Best');
    xlabel('Power Ratio');
    ylabel('empirical CDF');
    name = sprintf('results_mmse_ratio_%i',Rates(pl));
    saveas(gcf,sprintf('figures/%s.eps',name),'epsc2')
    saveas(gcf,sprintf('figures/%s.fig',name))
    
end