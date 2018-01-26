load('matlabPorownanieDMC.mat');

figure; plot((0:(length(h_out)-1))*Tp,h_out,'Color','black');
hold on;    plot((0:(length(hsp)-1))*Tp,hsp,'Color','b');
hold on;    plot((0:(length(h_out1)-1))*Tp,h_out1,'Color','red');
legend('h_o_u_t Analityczny','hsp','h_o_u_t Numeryczny','Location','northeast')
title('Wysokosc slupa cieczy')
xlabel('t[s]');
ylabel('h_o_u_t [cm]');

figure; plot((0:(length(Tout_v)-1))*Tp,Tout_v,'Color','black');
hold on;    plot((0:(length(Tvsp)-1))*Tp,Tvsp,'Color','b');
hold on;    plot((0:(length(Tout_v1)-1))*Tp,Tout_v1,'Color','red');
legend('T_o_u_t Analityczny','Tsp','T_o_u_t Numeryczny','Location','northeast')
title('Temperatura wyjsciowa')
xlabel('t[s]');
ylabel('T_o_u_t [C]');