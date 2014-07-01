P_tot_new = P_tot;
R_ac_new = R_ac;
R_ch_new = R_ch;
R_sum_new = R_sum;

%%

P_tot(:,6,:) = P_tot_new;
R_ac(:,6,:) = R_ac_new;
R_ch(:,:,6,:) = R_ch_new;
R_sum(:,6,:) = R_sum_new;

%% plots

figure
hold all
[channels, iterations] = size(P_now);
for ch = 1:channels
    plot(P_now(ch,:),'-','LineWidth',3);
end
title(sprintf('convergence of Powers for the numerical Gradient at %idB SNR',1));
xlabel('iteration');
ylabel('power');

figure
hold all
[channels, iterations] = size(P_now);
% for ch = 1:channels
%     plot(SINR_now(ch,:),'-','LineWidth',3);
% end
% SINRtgt = SINRtgt(1)*ones(1,iterations);
plot(sum_Rate,'-','LineWidth',3);%,'Color',[0 0 0]);
title(sprintf('convergence of the sum-rate for the numerical Gradient at %idB SNR',1));
xlabel('iteration');
ylabel('sum-rate');
