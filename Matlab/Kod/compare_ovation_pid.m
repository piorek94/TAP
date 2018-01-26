close all;
load('pid_decoupling_data.mat');
load pid_h_obrobione.csv

t_ov = pid_h_obrobione(:,1);
T_ov = pid_h_obrobione(:,3);
h_ov = pid_h_obrobione(:,2);
h_sp_ov = pid_h_obrobione(:,4);

model_step_time = 500;
ovation_step_time = 90;

pid_decoupling_data(1).Time = pid_decoupling_data(1).Time - model_step_time;
t_ov = t_ov - ovation_step_time;

fig = figure;
subplot(2,1,1);
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, (1)));
hold on;
plot(t_ov, h_ov);
hold on;
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, (4)));
xlim([-100 2500]);
title('Porownanie dzialania regulatorow PID - WYSOKOSC')
xlabel('t[s]');
ylabel('Wysokosc');
legend('Matlab', 'Ovation', 'h_{SP}', 'Location', 'southeast');
grid;

subplot(2,1,2);
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, (5)));
hold on;
plot(t_ov, T_ov);
hold on;
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, (8)));
xlim([-100 2500]);
title('Porownanie dzialania regulatorow PID - TEMPERATURA')
xlabel('t[s]');
ylabel('Temperatura');
legend('Matlab', 'Ovation', 'T_{SP}', 'Location', 'southeast');
grid;

print(fig, '../Dokumentacja/Obrazki/Etap3/pid_porownanie','-dpng','-r150');