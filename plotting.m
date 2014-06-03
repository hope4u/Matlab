%% Generate Plots for MIMO Sem Project Marc M�ller
close all
%% cdfs of the achievavle Rate for wf-Solutions
Typ = 2;
plots = {'none','wf','sp_iwf','sp_iwf_paper','numericalGrad_VBlast'};
names = {'none','waterfilling','iterative waterfilling','iterative waterfilling by jindal','numericalGradient'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1;1 1 0];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('water filling solutions at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_ac(snr,Typ,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
end

hold off
grid on
legend(leg);
end

%% cdfs for throughput Maximization
Typ = 1;
plots = {'none','numericalGrad','minmax','minmax_VBlast'};
names = {'none','numerical Gradient','max minimum Rate','max minimum Rate for VBlast'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('sum Rate Maximization at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if strcmp(plots(num),'minmax_VBlast')
        Typ = 2;
    end
    if Opt
        Vector = squeeze(R_sum(snr,Typ,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');   
        leg{end+1}=names{num};
    end
    if strcmp(plots(num),'minmax_VBlast')
        Typ = 1;
    end
end

hold off
grid on
legend(names);
end

%% cdfs for minimum Rates
Typ = 1;
plots = {'none','numericalGrad','minmax'};
names = {'none','numerical Gradient','max minimum Rate'};
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('minimum/maximum Rates at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_min(snr,Typ,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));

        Vector = squeeze(R_max(snr,Typ,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'--','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');
        leg{end+1}=sprintf('minRate of %s',names{num});
        leg{end+1}=sprintf('maxRate of %s',names{num});
    end
end

hold off
grid on
legend(leg);
end

%% cdfs for power Minimization
Typ = 1;
plots = {'powermin','fodorPrecoding2'};
names = {'numerical Gradient','Algorithm by Fodor'};
color = [1 0 0; 0 1 0; 0 0 1];

for snr=1:numel(SNR)
figure
hold on
title(sprintf('Power Minimization for Fixed Target SINR of 1dB at %i dB SNR',SNR(snr)));
leg = {};

for num=1:numel(plots)
    Opt = find(ismember(Optimizer,plots(num)));
    if Opt
        Vector = squeeze(R_sum(snr,Typ,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'-','LineWidth',3,'Color',color(num,:));

        Vector = squeeze(P_tot(snr,Typ,Opt,:));
        [f,x] = ecdf(Vector);
        plot(x,f,'--','LineWidth',3,'Color',color(num,:));
        set(gca, 'Yscale', 'log');
        leg{end+1}=names{num};
        leg{end+1}=sprintf('total Power');
    end
end

hold off
grid on
legend(leg);
end