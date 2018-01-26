close all;
load('tf_Fc_h_data.mat');
load('tf_Fc_T_data.mat');

model_step_time = 1500;
ovation_step_time = 1000;

load skok_Fc_obrobione.csv
t_ov = skok_Fc_obrobione(:,1);
T_ov = skok_Fc_obrobione(:,2);
h_ov = skok_Fc_obrobione(:,3);

tf_Fc_T_data(:,1) = tf_Fc_T_data(:,1) - model_step_time;
tf_Fc_h_data(:,1) = tf_Fc_h_data(:,1) - model_step_time;
t_ov = t_ov - ovation_step_time;

fig = figure;
plot(t_ov,T_ov);
hold on;
plot(tf_Fc_T_data(:,1), tf_Fc_T_data(:,2));
xlabel('t[s]');
ylabel('Temperatura');
title('Porownanie odpowiedzi skokowych modeli (F_C z 60 do 70)- TEMPERATURA');
legend('Ovation', 'Matlab - transmitancja');
xlim([0 750]);
grid;
print(fig, '../Dokumentacja/Obrazki/Etap3/skok_Fc_temperatura','-dpng','-r150');

fig = figure;
plot(t_ov,h_ov);
hold on;
plot(tf_Fc_h_data(:,1), tf_Fc_h_data(:,2));
title('Porownanie odpowiedzi skokowych modeli (F_C z 60 do 70)- WYSOKOSC');
xlabel('t[s]');
ylabel('Wysokosc');
legend('Ovation', 'Matlab - transmitancja');
xlim([-100 750]);
grid;
print(fig, '../Dokumentacja/Obrazki/Etap3/skok_Fc_wysokosc','-dpng','-r150');