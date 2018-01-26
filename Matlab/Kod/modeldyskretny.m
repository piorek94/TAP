%Projekt TAP - zbiornik z mieszadlem zadanie 3
%etap pierwszy podpunkt d - modele liniowe dyskretne - rownania stanu i
%transmitancja, sprawdzenie jakosci dyskretyzacji, porownanie modeli w
%przestrzeni stanu ciaglego i dyskretnego, porownanie odpowiedzi skokowych
%transmitancji ciaglej oraz dyskretnej
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
dFh=30;
dTh=30;
dFc=30;
dTc=30;
dFd=30;
dTd=30;
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
tauh0=230;
tau0=270;
%warunki poczatkowe
h0=18.65;
T0=33.16;

%--------------------------------------------------------------------------
%macierze przestrzeni stanu bez opoznien
Aps=[-alfa/(2*A*sqrt(h0)),0;-(Fh0*Th0+Fc0*Tc0+Fd0*Td0-(Fh0+Fc0+Fd0)*T0)/(A*h0^2),-(Fh0+Fc0+Fd0)/(A*h0)];
Bps=[1/A,0,1/A,0,1/A,0;(Th0-T0)/(A*h0),(Fh0)/(A*h0),(Tc0-T0)/(A*h0),(Fc0)/(A*h0),(Td0-T0)/(A*h0),(Fd0)/(A*h0)];
Cps=[1,0;0,1];
Dps=zeros(2,6);
%macierze przestrzeni stanu bez opoznien rozszerzona
Aps1=[-alfa/(2*A*sqrt(h0)),0,1;-(Fh0*Th0+Fc0*Tc0+Fd0*Td0-(Fh0+Fc0+Fd0)*T0)/(A*h0^2),-(Fh0+Fc0+Fd0)/(A*h0),0;0,0,0];
Bps1=[1/A,0,1/A,0,1/A,0;(Th0-T0)/(A*h0),(Fh0)/(A*h0),(Tc0-T0)/(A*h0),(Fc0)/(A*h0),(Td0-T0)/(A*h0),(Fd0)/(A*h0);0,0,0,0,0,0];
Cps1=[1,0,0;0,1,0;0,0,0];
Dps1=zeros(3,6);
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
sys_op=delayss(A0,B0,C0,D0,DelayT);
%model w przestrzeni stanu bez opoznien
sys=ss(Aps,Bps,Cps,Dps);
%model w przestrzeni stanu bez opoznien rozszerzony o stala
sys_1=ss(Aps1,Bps1,Cps1,Dps1);
%Mdl = ssm(Aps,Bps,Cps,Dps,'Mean0',[h0,T0]);

%--------------------------------------------------------------------------
%transmitancja z opoznieniami
G_op=tf(sys_op);
%transmitancja bez opoznien
G=tf(sys);
%transmitancja bez opoznien z modelu rozszerzonego
G_1=tf(sys_1);
%--------------------------------------------------------------------------

%czas probkowania oraz metoda dyskretyzacji
Tp=10;
metoda='tustin';


%modele dyskretne

%modele w postaci transmitancji dyskretnej
%model dyskretny liniowy w postaci transmitancji z opoznieniami
Gz_op=c2d(G_op,Tp,metoda);
%model dyskretny liniowy w postaci transmitancji bez opoznien
Gz=c2d(G,Tp,metoda);
%model dyskretny liniowy w postaci transmitancji bez opoznien rozszerzony z
%modelu rozszerzonego
Gz_1=c2d(G_1,Tp,metoda);

%model w przestrzeni stanu dyskretny
%model dyskretny liniowy w postaci rownan stanu z opoznieniami
sysz_op=c2d(sys_op,Tp,metoda);
%model dyskretny liniowy w postaci rownan stanu bez opoznien
sysz=c2d(sys,Tp,metoda);
%model dyskretny liniowy w postaci rownan stanu rozszerzony bez opoznien
sysz_1=c2d(sys_1,Tp,metoda);
[Az,Bz,Cz,Dz]=ssdata(sysz_1);

%--------------------------------------------------------------------------
%porownanie okresu probkowania
for i=1:6
    for j=1:2
        Tp1=50; % pocz¹tkowa wartoœæ okresu próbkowania (wyraŸnie zbyt du¿a)
        figure
        sysd=c2d(sys_op,Tp1,metoda);
        step(sys_op(j,i)); 
        hold on; 
        step(sysd(j,i))

        Tp1=Tp1/5;%10
        sysd1=c2d(sys_op,Tp1,metoda);
        step(sysd1(j,i));

        Tp1=Tp1/2;%5
        sysd4=c2d(sys_op,Tp1,metoda);
        step(sysd4(j,i));

        Tp1=Tp1/5;%1
        sysd5=c2d(sys_op,Tp1,metoda);
        step(sysd5(j,i));
        %print(num2str(j*10+i),'-dpng','-r600');
    end
end

%--------------------------------------------------------------------------
%porownanie modelu w przestrzeni stanu ciaglego i dyskretnego
for i=1:6
figure

switch i
    case 1
        a1=1;
        a2=0;
        a3=0;
        a4=0;
        a5=0;
        a6=0;
    case 2
        a1=0;
        a2=1;
        a3=0;
        a4=0;
        a5=0;
        a6=0;
    case 3
        a1=0;
        a2=0;
        a3=1;
        a4=0;
        a5=0;
        a6=0;
    case 4
        a1=0;
        a2=0;
        a3=0;
        a4=1;
        a5=0;
        a6=0;
    case 5
        a1=0;
        a2=0;
        a3=0;
        a4=0;
        a5=1;
        a6=0;
    case 6
        a1=0;
        a2=0;
        a3=0;
        a4=0;
        a5=0;
        a6=1;
end

%skoki
deltaFh=dFh*a1;
deltaTh=dTh*a2;
deltaFc=dFc*a3;
deltaTc=dTc*a4;
deltaFd=dFd*a5;
deltaTd=dTd*a6;

%wywolanie modelu na podstawaie modelu ciaglego w przestrzeni stanu
sim('prz_st',Tsym)

%wywolanie modelu na podstawaie modelu dyskretnego w przestrzeni stanu
sim('przestrzen_stanu_dys',Tsym)
    
%charakterystyka wysokosci    

subplot(1,2,1);
stairs(ps_dys_h.Time,ps_dys_h.Data,'Color','b','LineWidth',0.25)
hold on    
plot(prz_st_h.Time,prz_st_h.Data,'LineStyle',':','Color','r','LineWidth',0.75)

legend('przestrzen stanu dyskretny','przestrzen stanu ciagly','Location','southeast')
title('wysokosc h')
xlabel('t[s]');
ylabel('h[cm]');

%charakterystyka temperatury
subplot(1,2,2);
stairs(ps_dys_tout.Time,ps_dys_tout.Data,'Color','b','LineWidth',0.25)
hold on
plot(prz_st_tout.Time,prz_st_tout.Data,'LineStyle',':','Color','r','LineWidth',0.75)

legend('przestrzen stanu dyskretny','przestrzen stanu ciagly','Location','southeast')
title('temperatura Tout')
xlabel('t[s]');
ylabel('T_o_u_t [C]');

%print(num2str(i),'-dpng','-r600');
end

%{
%porownanie modelu w przestrzeni stanu ciaglego i dyskretnego
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

    %wywolanie modelu na podstawaie modelu ciaglego w przestrzeni stanu
    sim('prz_st',Tsym)

    %wywolanie modelu na podstawaie modelu dyskretnego w przestrzeni stanu
    sim('przestrzen_stanu_dys',Tsym)
    
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
    stairs(ps_dys_h.Time,ps_dys_h.Data,'Color',col,'LineWidth',0.25)
    hold on    
    plot(prz_st_h.Time,prz_st_h.Data,'Color',col,'LineStyle',':','Marker','o','MarkerSize',1,'LineWidth',0.25)
    hold on
end
%--------------------------------------------------------------------------

legend('przestrzen stanu ciagly','przestrzen stanu dyskretny','Location','southeast')
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

    %wywolanie modelu na podstawaie modelu ciaglego w przestrzeni stanu
    sim('prz_st',Tsym)

    %wywolanie modelu na podstawaie modelu dyskretnego w przestrzeni stanu
    sim('przestrzen_stanu_dys',Tsym)
    
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
    stairs(ps_dys_tout.Time,ps_dys_tout.Data,'Color',col,'LineWidth',0.25)
    hold on
    plot(prz_st_tout.Time,prz_st_tout.Data,'Color',col,'LineStyle',':','Marker','o','MarkerSize',1,'LineWidth',0.25)
    hold on
end

legend('przestrzen stanu ciagly','przestrzen stanu dyskretny','Location','southeast')
title('temperatura Tout')
xlabel('t[s]');
ylabel('T_o_u_t [C]');
print(nazwa,'-dpng','-r600');
%}
%porownanie odpowiedzi skokowych transmitancji ciaglej i dyskretnej
figure
step(G_op);hold on; step(Gz_op);
%print('odp_skok_tran_cia_dys','-dpng','-r600');

