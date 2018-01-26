close all;
load('tf_Fh_h_data.mat');
load('tf_Fh_T_data.mat');

model_step_time = 1500;
ovation_step_time = 500;

load skok_Fh_obrobione.csv
t_ov = skok_Fh_obrobione(:,1);
T_ov = skok_Fh_obrobione(:,2);
h_ov = skok_Fh_obrobione(:,3);

tf_Fh_T_data(:,1) = tf_Fh_T_data(:,1) - model_step_time;
tf_Fh_h_data(:,1) = tf_Fh_h_data(:,1) - model_step_time;
t_ov = t_ov - ovation_step_time;

fig = figure;
plot(t_ov,T_ov);
hold on;
plot(tf_Fh_T_data(:,1), tf_Fh_T_data(:,2));
xlabel('t[s]');
ylabel('Temperatura');
title('Porownanie odpowiedzi skokowych modeli (F_H z 20 do 30)- TEMPERATURA');
legend('Ovation', 'Matlab - transmitancja');
xlim([0 1500]);
grid;
print(fig, '../Dokumentacja/Obrazki/Etap3/skok_Fh_temperatura','-dpng','-r150');

fig = figure;
plot(t_ov,h_ov);
hold on;
plot(tf_Fh_h_data(:,1), tf_Fh_h_data(:,2));
title('Porownanie odpowiedzi skokowych modeli (F_H z 20 do 30)- WYSOKOSC');
xlabel('t[s]');
ylabel('Wysokosc');
legend('Ovation', 'Matlab - transmitancja');
xlim([0 1500]);
grid;
print(fig, '../Dokumentacja/Obrazki/Etap3/skok_Fh_wysokosc','-dpng','-r150');