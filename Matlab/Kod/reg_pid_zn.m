close all;
decouplingEnabled = 0;
ki_h = 0;
ki_T = 0;
Tsym = 12000;

%dobor nastaw regulatora wysokosci
deltah = 10;
deltaT = 0;
kp_h = 5.24;
kp_T = 0;
sim('pid_a');
figh = figure;
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 1));
%hold on;
%plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, 1));
title('Regulacja wysokosci');
%legend('bez odsprzegania', 'z odsprzeganiem');
grid;

%dobor nastaw regulatora temperatury
deltah = 0;
deltaT = 10;
kp_h = 0;
kp_T = 7.9;
sim('pid_a');

figT = figure;
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 5));
%hold on;
%plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, 5));
title('Regulacja temperatury');
%legend('bez odsprzegania', 'z odsprzeganiem');
grid;
