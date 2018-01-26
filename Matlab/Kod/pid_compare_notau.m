% h and T control comparison on one figure
fighT_comparison_notau = figure;
subplot(2,1,1);
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 1));
hold on;
plot(pid_notau_data(1).Time, pid_notau_data(1).Data(:, 1));
hold on;
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 4), '--');
title('Regulacja wysokosci, porownanie');
legend('z opoznieniem tau=270s', 'z opoznieniem tau=0s');
grid;

subplot(2,1,2);
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 5));
hold on;
plot(pid_notau_data(1).Time, pid_notau_data(1).Data(:, 5));
hold on;
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, 8), '--');
title('Regulacja temperatury, porownanie');
legend('z opoznieniem tau=270s', 'z opoznieniem tau=0s');
grid;
print(fighT_comparison_notau, '../Dokumentacja/Obrazki/Etap2/pid_hT_comparison_notau', '-dpng', '-r150');