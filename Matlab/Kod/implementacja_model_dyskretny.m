%Projekt TAP - zbiornik z mieszadlem zadanie 3
%etap pierwszy podpunkd d-impelementacja modelu dyskretnego -rownania stanu
%porównanie modelu dyskretnego zaimplementowanego i ciaglego
close all
clear
clc
%stale
alfa=22;
A=500;
%czas symulowany
Tsym=4000;
Tskok=1500;
Tciagly=1;
%--------------------------------------------------------------------------
%nazwy rysunkow
nazwa='porownanie-Fh';
%skoki
dFh=10;
dTh=0;
dFc=0;
dTc=0;
dFd=0;
dTd=0;
%sygnaly
Fh=20;
Th=75;
Fc=60;
Tc=17;
Fd=15;
Td=42;
%opoznienia w modelach
tauh=230;
tau=270;
h=18.65;
T=33.16;
%--------------------------------------------------------------------------
%Punkt pracy/linearyzacji
Tc0=17;
Th0=75;
Td0=42;
Fc0=60;
Fh0=20;
Fd0=15;
%warunki poczatkowe
h0=18.65;
T0=33.16;
%macierze przestrzeni stanu bez opoznien rozszerzona
Aps1=[-alfa/(2*A*sqrt(h0)),0,1;-(Fh0*Th0+Fc0*Tc0+Fd0*Td0-(Fh0+Fc0+Fd0)*T0)/(A*h0^2),-(Fh0+Fc0+Fd0)/(A*h0),0;0,0,0];
Bps1=[1/A,0,1/A,0,1/A,0;(Th0-T0)/(A*h0),(Fh0)/(A*h0),(Tc0-T0)/(A*h0),(Fc0)/(A*h0),(Td0-T0)/(A*h0),(Fd0)/(A*h0);0,0,0,0,0,0];
Cps1=[1,0,0;0,1,0;0,0,0];
Dps1=zeros(3,6);
%model w przestrzeni stanu bez opoznien rozszerzony o stala
sys_1=ss(Aps1,Bps1,Cps1,Dps1);
%czas probkowania oraz metoda dyskretyzacji
Tp=10;
metoda='tustin';
%model dyskretny liniowy w postaci rownan stanu rozszerzony bez opoznien
sysz_1=c2d(sys_1,Tp,metoda);
[Az,Bz,Cz,Dz]=ssdata(sysz_1);

Fh_v=zeros(1,Tsym/Tp);
Fh_v=Fh_v+Fh;
Fh_v(1,(Tskok+tauh)/Tp +1:end)=Fh_v(1,(Tskok+tauh)/Tp +1:end)+dFh;

Th_v=zeros(1,Tsym/Tp);
Th_v=Th_v+Th;
Th_v(1,(Tskok)/Tp +1:end)=Th_v(1,(Tskok)/Tp +1:end)+dTh;

Fd_v=zeros(1,Tsym/Tp);
Fd_v=Fd_v+Fd;
Fd_v(1,(Tskok)/Tp +1:end)=Fd_v(1,(Tskok)/Tp +1:end)+dFd;

Td_v=zeros(1,Tsym/Tp);
Td_v=Td_v+Td;
Td_v(1,(Tskok)/Tp +1:end)=Td_v(1,(Tskok)/Tp +1:end)+dTd;

Fc_v=zeros(1,Tsym/Tp);
Fc_v=Fc_v+Fc;
Fc_v(1,(Tskok)/Tp +1:end)=Fc_v(1,(Tskok)/Tp +1:end)+dFc;

Tc_v=zeros(1,Tsym/Tp);
Tc_v=Tc_v+Tc;
Tc_v(1,(Tskok)/Tp +1:end)=Tc_v(1,(Tskok)/Tp +1:end)+dTc;

T_v=zeros(1,Tsym/Tp);
T_v=T_v+T0;

h_v=zeros(1,Tsym/Tp);
h_v=h_v+h0;
h_out = h_v;

Tout_v=zeros(1,Tsym/Tp);
Tout_v=Tout_v+T0;


for i=2:Tsym/Tp
    h_v(i)=Az(1,1)*h_v(i-1)+Bz(1,1)*Fh_v(i-1)+Bz(1,3)*Fc_v(i-1)+Bz(1,5)*Fd_v(i-1)+Az(1,3)*(-alfa*sqrt(h0)/(A*2));
    h_out(i)=Cz(1,1)*h_v(i)+Dz(1,1)*Fh_v(i)+Dz(1,3)*Fc_v(i)+Dz(1,5)*Fd_v(i)+Cz(1,3)*(-alfa*sqrt(h0)/(A*2));
    
    T_v(i)=Az(2,1)*h_v(i-1)+Az(2,2)*T_v(i-1)+Bz(2,1)*Fh_v(i-1)+Bz(2,2)*Th_v(i-1)+Bz(2,3)*Fc_v(i-1)+Bz(2,4)*Tc_v(i-1)+Bz(2,5)*Fd_v(i-1)+Bz(2,6)*Td_v(i-1) + Az(2,3)*(-alfa*sqrt(h0)/(A*2));
    if(i>tau/Tp)
        Tout_v(i)=Cz(2,1)*h_v(i-tau/Tp)+Cz(2,2)*T_v(i-tau/Tp)+Dz(2,1)*Fh_v(i-tau/Tp)+Dz(2,2)*Th_v(i-tau/Tp)+Dz(2,3)*Fc_v(i-tau/Tp)+Dz(2,4)*Tc_v(i-tau/Tp)+Dz(2,5)*Fd_v(i-tau/Tp)+Dz(2,6)*Td_v(i-tau/Tp)+Cz(2,3)*(-alfa*sqrt(h0)/(A*2));
    end
    czas(i)=(i-1)*Tp;
end
deltaFh=dFh;
deltaTh=dTh;
deltaFc=dFc;
deltaTc=dTc;
deltaFd=dFd;
deltaTd=dTd;
%wywolanie modelu ciaglego w przestrzeni stanu
sim('prz_st',Tsym)

%wysokosc
subplot(1,2,1);
stairs(czas,h_out,'Color','b','LineWidth',0.25)
hold on
plot(prz_st_h.Time,prz_st_h.Data,'LineStyle',':','Color','r','LineWidth',0.75)

legend('przestrzen stanu dyskretny','przestrzen stanu ciagly','Location','southeast')
title('wysokosc h')
xlabel('t[s]');
ylabel('h[cm]');

%temperatura
subplot(1,2,2);
stairs(czas,Tout_v,'Color','b','LineWidth',0.25)
hold on
plot(prz_st_tout.Time,prz_st_tout.Data,'LineStyle',':','Color','r','LineWidth',0.75)

legend('przestrzen stanu dyskretny','przestrzen stanu ciagly','Location','southeast')
title('temperatura Tout')
xlabel('t[s]');
ylabel('T_o_u_t [C]');

%print('por_ciag_dysk_imple','-dpng','-r600');