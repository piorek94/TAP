close all;
load('pid_decoupling_data.mat');
load('dmc_data.mat');

fig = figure;
subplot(2,1,1);
p1 = plot(pid_decoupling_data(1).Time, [pid_decoupling_data(1).Data(:, (1)) pid_decoupling_data(1).Data(:, (4))]);
p1(1).LineWidth = 2;
p1(2).LineWidth = 2;
hold on;
p2 = plot((0:(length(h_out)-1))*Tp,h_out);
p2.LineWidth = 2;
legend('h_{PID}', 'h_{sp}', 'h_{DMC}');
lowerTitle = 'i analitycznego DMC, skok zaklocenia $F_d$ w t=12000s';
title({'\makebox[12cm][c]{Porownanie dzialania regulatorow wysokosci: PID z odsprzeganiem}',strcat('\makebox[12cm][c]{', lowerTitle, '}')},'Interpreter','latex');
xlabel('Czas [s]');
grid;
subplot(2,1,2);
p3 = plot(pid_decoupling_data(1).Time, [pid_decoupling_data(1).Data(:, (5)) pid_decoupling_data(1).Data(:, (8))]);
p3(1).LineWidth = 2;
p3(2).LineWidth = 2;
hold on;
p4 = plot((0:(length(Tout_v)-1))*Tp,Tout_v); %,'Color','black'
p4.LineWidth = 2;
legend('T_{PID}', 'T_{sp}', 'T_{DMC}');
lowerTitle = 'i analitycznego DMC, skok zaklocenia $F_d$ w t=12000s';
title({'\makebox[12cm][c]{Porownanie dzialania regulatorow temperatury: PID z odsprzeganiem}',strcat('\makebox[12cm][c]{', lowerTitle, '}')},'Interpreter','latex')
xlabel('Czas [s]');
grid;

print(fig, '../Dokumentacja/Obrazki/Etap2/comparison_pid_dmc', '-dpng', '-r150');