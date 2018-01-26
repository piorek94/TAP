%TAP3_1; % zaladowanie transmitancji, punktow pracy itd.
close all;
Tsym = 16000;
Tp_pid = 1;
deltah = 5; % skok wartosci zadanej h
deltaT = 5; % skok wartosci zadanej T
deltaFd = 5; % skok zaklocenia
h_step_up_time = 500;
h_step_down_time = 8000;
T_step_up_time = 4000;
T_step_down_time = 8000;
disturbanceTime = 12000;

% wlaczanie transmitancji skrosnych
enable_cross_tf = 1;

%nastawy regulatorow - bez odsprzegania
kp_h_normal = 5;
Ti_h_normal = 50; %1 / 0.01;
kd_h_normal = 4;
taud_h_normal = 2;

kp_T_normal = 1;
Ti_T_normal = 400;
kd_T_normal = 0;
taud_T_normal = 0.5;

%nastawy regulatorow - z odsprzeganiem
kp_h_decoupling = 5;
Ti_h_decoupling = 50;
kd_h_decoupling = 4;
taud_h_decoupling = 2;

kp_T_decoupling = 1;
Ti_T_decoupling = 400;
kd_T_decoupling = 0;
taud_T_decoupling = 0.5;

%zakresy wartosci sterowanych, wzgledem punktu pracy Fc i Fh
enable_saturation = 1;
Fc_min = 0 - Fc;
Fc_max = 100 - Fc;
Fh_min = 0 - Fh;
Fh_max = 100 - Fh;
if(enable_saturation == 0)
   Fc_max = 999999;
   Fh_max = 999999;
   Fc_min = -Fc_max;
   Fh_min = -Fh_max;
end


decouplingEnabled = 1;
pid_ovation_h = struct('Tp', 1, 'Kp', kp_h_decoupling, 'Ti', Ti_h_decoupling, 'Kd', kd_h_decoupling, 'taud', taud_h_decoupling);
pid_ovation_T = struct('Tp', 1, 'Kp', kp_T_decoupling, 'Ti', Ti_T_decoupling, 'Kd', kd_T_decoupling, 'taud', taud_T_decoupling);
D_12 = -G(1,1) / G(1,3);
%D_12 = tf(-1, 1);
D_12.IOdelay = G1(1,1).IOdelay - G1(1,3).IOdelay;

D_21 = -G(2,3) / G(2,1) * 0;
%D_21 = tf(0.0776, 1);

sim('regulacja_pid');
pid_decoupling_data = pid_a_data;

decouplingEnabled = 0;
pid_ovation_h = struct('Tp', 1, 'Kp', kp_h_normal, 'Ti', Ti_h_normal, 'Kd', kd_h_normal, 'taud', taud_h_normal);
pid_ovation_T = struct('Tp', 1, 'Kp', kp_T_normal, 'Ti', Ti_T_normal, 'Kd', kd_T_normal, 'taud', taud_T_normal);
sim('regulacja_pid');


% all h without decoupling variables
% figh = figure;
% % plot(pid_a_data(1).Time, pid_a_data(1).Data(:, [1:4]));
% % legend('h', 'F_C', 'e_h', 'h_{sp}');
% plot(pid_a_data(1).Time, [pid_a_data(1).Data(:, (1)) pid_a_data(1).Data(:, (4))]);
% legend('h', 'h_{sp}');
% title('Regulacja wysokosci bez transmitancji skrosnych');
% xlabel('Czas [s]');
% grid;
%print(figh, '../Dokumentacja/Obrazki/Etap2/pid_h_nodec_nocross', '-dpng', '-r150');

% all T without decoupling variables
% figT = figure;
% % plot(pid_a_data(1).Time, pid_a_data(1).Data(:, [5:8]));
% % legend('T', 'F_H', 'e_T', 'T_{sp}');
% plot(pid_a_data(1).Time, [pid_a_data(1).Data(:, (5)) pid_a_data(1).Data(:, (8))]);
% legend('T', 'T_{sp}');
% title('Regulacja temperatury bez transmitancji skrosnych');
% xlabel('Czas [s]');
% grid;
%print(figT, '../Dokumentacja/Obrazki/Etap2/pid_T_nodec_nocross', '-dpng', '-r150');

% h and T without decoupling on one figure
fighT = figure;
subplot(2,1,1);
plot(pid_a_data(1).Time, [pid_a_data(1).Data(:, (1)) pid_a_data(1).Data(:, (4))]);
legend('h', 'h_{sp}');
title('Regulacja wysokosci bez osprzegania, skok zaklocenia F_d w t=12000s');
xlabel('Czas [s]');
grid;
subplot(2,1,2);
plot(pid_a_data(1).Time, [pid_a_data(1).Data(:, (5)) pid_a_data(1).Data(:, (8))]);
legend('T', 'T_{sp}');
title('Regulacja temperatury bez odsprzegania, skok zaklocenia F_d w t=12000s');
xlabel('Czas [s]');
grid;
%print(fighT, '../Dokumentacja/Obrazki/Etap2/pid_hT_nodec', '-dpng', '-r150');

% h and T with decoupling on one figure
fighT_dec = figure;
subplot(2,1,1);
plot(pid_decoupling_data(1).Time, [pid_decoupling_data(1).Data(:, (1)) pid_decoupling_data(1).Data(:, (4))]);
legend('h', 'h_{sp}');
title('Regulacja wysokosci z odsprzeganiem, skok zaklocenia F_d w t=12000s');
xlabel('Czas [s]');
grid;
subplot(2,1,2);
plot(pid_decoupling_data(1).Time, [pid_decoupling_data(1).Data(:, (5)) pid_decoupling_data(1).Data(:, (8))]);
legend('T', 'T_{sp}');
title('Regulacja temperatury z odsprzeganiem, skok zaklocenia F_d w t=12000s');
xlabel('Czas [s]');
grid;
%print(fighT_dec, '../Dokumentacja/Obrazki/Etap2/pid_hT_dec', '-dpng', '-r150');



%all h with decoupling variables
% figh_dec = figure;
% plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, [1:4]));
% legend('h', 'F_H', 'e_h', 'h_{sp}');
% % plot(pid_decoupling_data(1).Time, [pid_decoupling_data(1).Data(:, (1)) pid_decoupling_data(1).Data(:, (4))]);
% % legend('h', 'h_{sp}');
% title('Regulacja wysokosci z odsprzeganiem');
% xlabel('Czas [s]');
% grid;

% all T with decoupling variables
% figT_dec = figure;
% plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, [5:8]));
% legend('T', 'F_C', 'e_T', 'T_{sp}');
% % plot(pid_decoupling_data(1).Time, [pid_decoupling_data(1).Data(:, (5)) pid_decoupling_data(1).Data(:, (8))]);
% % legend('T', 'T_{sp}');
% title('Regulacja temperatury z odsprzeganiem');
% xlabel('Czas [s]');
% grid;

% h and T control comparison on one figure
fighT_comparison = figure;
subplot(2,1,1);
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 1));
hold on;
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, 1));
hold on;
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 4), '--');
title('Regulacja wysokosci, porownanie');
legend('bez odsprzegania', 'z odsprzeganiem');
grid;

subplot(2,1,2);
plot(pid_a_data(1).Time, pid_a_data(1).Data(:, 5));
hold on;
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, 5));
hold on;
plot(pid_decoupling_data(1).Time, pid_decoupling_data(1).Data(:, 8), '--');
title('Regulacja temperatury, porownanie');
legend('bez odsprzegania', 'z odsprzeganiem');
grid;
%print(fighT_comparison, '../Dokumentacja/Obrazki/Etap2/pid_hT_comparison', '-dpng', '-r150');