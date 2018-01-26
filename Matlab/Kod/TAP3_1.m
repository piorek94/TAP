%Projekt TAP - zbiornik z mieszadlem zadanie 3
%etap pierwszy podpunkt a, porownanie modeli z rzeczywistymi warunakmi
%poczatkowymi
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
%skoki
deltaFh=10;
deltaTh=0;
deltaFc=0;
deltaTc=0;
deltaFd=0;
deltaTd=0;
%sygnaly
Fh=20;
Th=75;
Fc=70;
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
tauh0=230;
tau0=270;
%warunki poczatkowe
h0=18.65;
T0=33.16;
%rzeczywiste warunki pocz�tkowe
h0=((Fh0+Fc0+Fd0)/(alfa))^2;
T0=(Fh0*Th0+Fc0*Tc0+Fd0*Td0)/(Fh0+Fc0+Fd0);

%--------------------------------------------------------------------------
%macierze przestrzeni stanu bez opoznien
Aps=[-alfa/(2*A*sqrt(h0)),0;-(Fh0*Th0+Fc0*Tc0+Fd0*Td0-(Fh0+Fc0+Fd0)*T0)/(A*h0^2),-(Fh0+Fc0+Fd0)/(A*h0)];
Bps=[1/A,0,1/A,0,1/A,0;(Th0-T0)/(A*h0),(Fh0)/(A*h0),(Tc0-T0)/(A*h0),(Fc0)/(A*h0),(Td0-T0)/(A*h0),(Fd0)/(A*h0)];
Cps=[1,0;0,1];
Dps=zeros(2,6);

%macierze do modelu z opoznieniami
A0=[-alfa/(2*A*sqrt(h0)),0;-(Fh0*Th0+Fc0*Tc0+Fd0*Td0-(Fh0+Fc0+Fd0)*T0)/(A*h0^2),-(Fh0+Fc0+Fd0)/(A*h0)];
A1=zeros(2,2);
A2=zeros(2,2);
D0=zeros(2,6);
D1=zeros(2,6);
D2=zeros(2,6);
C0=[1,0;0,0];
C1=zeros(2,2);
C2=[0,0;0,1];
B0=[0,0,1/A,0,1/A,0;0,(Fh0)/(A*h0),(Tc0-T0)/(A*h0),(Fc0)/(A*h0),(Td0-T0)/(A*h0),(Fd0)/(A*h0)];
B1=[1/A,0,0,0,0,0;(Th0-T0)/(A*h0),0,0,0,0,0];
B2=zeros(2,6);
%kolejne opoznienia
DelayT(1)=struct('delay',tauh,'a',A1,'b',B1,'c',C1,'d',D1);
DelayT(2)=struct('delay',tau,'a',A2,'b',B2,'c',C2,'d',D2);

%--------------------------------------------------------------------------
%model w przestrzeni stanu z opoznieniami
sys1=delayss(A0,B0,C0,D0,DelayT);

%model w przestrzeni stanu bez opoznien
sys=ss(Aps,Bps,Cps,Dps);
%Mdl = ssm(Aps,Bps,Cps,Dps,'Mean0',[h0,T0]);

%--------------------------------------------------------------------------
%transmitancja z opoznieniami
G1=tf(sys1);

%transmitancja bez opoznien
G=tf(sys);

%--------------------------------------------------------------------------
%wywo�anie modelu na podstawie rownan rozniczkowych nieliniowych
sim('roz_nieliniowe',Tsym)

%wywolanie modelu na podstawaie modelu w przestrzeni stanu
sim('prz_st',Tsym)

%wywolanie modelu na podstawaie transmitancji-brak warunkow poczatkowych i
%stalej
sim('tran_opoz',Tsym)

%rozwiazanie rownan rozniczkowych nieliniowych
[Tnl,Ynl] = ode45(@(t,y) sym1(t,y,Fh,Th,Fc,Tc,Fd,Td,alfa,A),[0 Tsym],[h0 T0]);

%rozwiazanie rownan rozniczkowych liniowych
[Tl,Yl] = ode45(@(t,y) sym2( t,y,Fh,Th,Fc,Tc,Fd,Td,alfa,A,Tc0,Th0,Td0,Fc0,Fh0,Fd0,h0,T0 ),[0 Tsym],[h0 T0]);


%--------------------------------------------------------------------------
%charakterystyka wysokosci

%z rownaniami rozniczkowymi
%plot(Tnl,Ynl(:,1),'-',Tl,Yl(:,1),'-.')%,ro_nl_h.Time,ro_nl_h.Data,ro_lin_h.Time,ro_lin_h.Data,'--',ps_h.Time,ps_h.Data,':',tran_h.Time,tran_h.Data)
%legend('nieliniowy','liniowy','model nieliniowy','model liniowy','przestrzen stanu','transmitancja','Location','southeast')

%bez rownan rozniczkowych
plot(ro_nl_h.Time,ro_nl_h.Data,prz_st_h.Time,prz_st_h.Data,':',tran_h.Time,tran_h.Data,'--')
legend('model nieliniowy','przestrzen stanu','transmitancja','Location','southeast')
title('wysokosc h')
xlabel('t[s]');
ylabel('h[cm]');
print('../Dokumentacja/Obrazki/wysokosc3','-dpng','-r600');


%--------------------------------------------------------------------------
%charakterystyka temperatury
figure
%z rownaniami rozniczkowymi
%plot(Tnl,Ynl(:,2),'-',Tl,Yl(:,2),'-.',ro_nl_tout.Time,ro_nl_tout.Data,ro_lin_tout.Time,ro_lin_tout.Data,'--',ps_tout.Time,ps_tout.Data,':',tran_tout.Time,tran_tout.Data)
%legend('nieliniowy','liniowy','model nieliniowy','model liniowy','przestrzen stanu','transmitancja','Location','southeast')

%bez rownan rozniczkowych
plot(ro_nl_tout.Time,ro_nl_tout.Data,prz_st_tout.Time,prz_st_tout.Data,':',tran_tout.Time,tran_tout.Data,'--')
legend('model nieliniowy','przestrzen stanu','transmitancja','Location','southeast')
title('temperatura Tout')
xlabel('t[s]');
ylabel('T_o_u_t [C]');
print('../Dokumentacja/Obrazki/temperatura3','-dpng','-r600');


%--------------------------------------------------------------------------
%model dyskretny liniowy w postaci transmitancji
Gz=c2d(G1,0.5,'zoh');
%model dyskretny liniowy w postaci rownan stanu
sysz=c2d(sys1,0.5,'zoh');