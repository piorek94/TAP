%Projekt TAP - zbiornik z mieszadlem zadanie 3
%etap pierwszy podpunkt b - porownanie  model liniowy ( przestrzen stanu i
%transmitancja ) z modelem nieliniowym dla roznych skokow wejsc
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
nazwa='wys_tem_';
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
%rzeczywiste warunki pocz¹tkowe
%h0=((Fh0+Fc0+Fd0)/(alfa))^2;
%T0=(Fh0*Th0+Fc0*Tc0+Fd0*Td0)/(Fh0+Fc0+Fd0);
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
%transmitancja z opoznieniami
G1=tf(sys1);
%--------------------------------------------------------------------------

figure
subplot(1,2,1);
%charakterystyka wysokosci
for i=0:3
%skoki
deltaFh=dFh*i;
deltaTh=dTh*i;
deltaFc=dFc*i;
deltaTc=dTc*i;
deltaFd=dFd*i;
deltaTd=dTd*i;
    %wywo³anie modelu na podstawie rownan rozniczkowych nieliniowych
    sim('roz_nieliniowe',Tsym)

    %wywolanie modelu na podstawaie modelu w przestrzeni stanu
    sim('prz_st',Tsym)

    %wywolanie modelu na podstawaie transmitancji-brak warunkow poczatkowych i
    %stalej
    sim('tran_opoz',Tsym)
    
    %charakterystyka wysokosci
    switch i
        case 0
    col='k';
        case 1
    col='r';   
        case 2
    col='b';
        case 3
    col='g';        
    end
    
    plot(ro_nl_h.Time,ro_nl_h.Data,'Color',col,'LineWidth',0.25)
    hold on
    plot(prz_st_h.Time,prz_st_h.Data,'Color',col,'LineStyle',':','Marker','o','MarkerSize',1,'LineWidth',0.25)
    hold on
    plot(tran_h.Time,tran_h.Data,'Color',col,'LineStyle','--','LineWidth',0.25)
    hold on    
end
%--------------------------------------------------------------------------

legend('model nieliniowy','przestrzen stanu','transmitancja','Location','southeast')
title('wysokosc h')
xlabel('t[s]');
ylabel('h[cm]');


%--------------------------------------------------------------------------
%charakterystyka temperatury
subplot(1,2,2);
for i=0:3
%skoki
deltaFh=dFh*i;
deltaTh=dTh*i;
deltaFc=dFc*i;
deltaTc=dTc*i;
deltaFd=dFd*i;
deltaTd=dTd*i;
    %wywo³anie modelu na podstawie rownan rozniczkowych nieliniowych
    sim('roz_nieliniowe',Tsym)

    %wywolanie modelu na podstawaie modelu w przestrzeni stanu
    sim('prz_st',Tsym)

    %wywolanie modelu na podstawaie transmitancji-brak warunkow poczatkowych i
    %stalej
    sim('tran_opoz',Tsym)
    
    %charakterystyka wysokosci
    switch i
        case 0
    col='k';
        case 1
    col='r';   
        case 2
    col='b';
        case 3
    col='g';        
    end
    plot(ro_nl_tout.Time,ro_nl_tout.Data,'Color',col,'LineWidth',0.25)
    hold on
    plot(prz_st_tout.Time,prz_st_tout.Data,'Color',col,'LineStyle',':','Marker','o','MarkerSize',1,'LineWidth',0.25)
    hold on
    plot(tran_tout.Time,tran_tout.Data,'Color',col,'LineStyle','--','LineWidth',0.25)
    hold on    
end

legend('model nieliniowy','przestrzen stanu','transmitancja','Location','southeast')
title('temperatura Tout')
xlabel('t[s]');
ylabel('T_o_u_t [C]');
print(nazwa,'-dpng','-r600');