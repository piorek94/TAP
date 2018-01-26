clear
clc
options = optimoptions('quadprog','Algorithm','interior-point-convex');
%%% Parametry poczatkowe modelu
alfa=22;
A=500;
%czas symulowany
Tsym=16000;
Tskok=11000;
Tskokzaklocenie=14000;
hSkokUp = 500;
hSkokDown = 6000;
TvSkokUp = 3000;
TvSkokDown = 9000;
Tciagly=1;
%skoki
dFh=0;
dTh=0;
dFc=0;
dTc=0;
dFd=10;
dTd=10;
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
%rzeczywiste warunki poczatkowe
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

h0=((Fh0+Fc0+Fd0)/(alfa))^2;
T0=(Fh0*Th0+Fc0*Tc0+Fd0*Td0)/(Fh0+Fc0+Fd0);

hsp = zeros(1,Tsym/Tp);
hsp = hsp+h0;
hsp(1,hSkokUp/Tp:hSkokDown/Tp)=h0+10;

Tvsp = zeros(1,Tsym/Tp);
Tvsp = Tvsp+T0;
Tvsp(1,TvSkokUp/Tp:TvSkokDown/Tp)=T0+10;
FH = zeros(1,Tsym/Tp);
FH=FH+Fh;

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
Td_v(1,(Tskokzaklocenie)/Tp +1:end)=Td_v(1,(Tskokzaklocenie)/Tp +1:end)+dTd;

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

%%% Parametry poczatkowe dla DMC
Tp = 10;  % czas probkowania
D = 160; % horyzont Dynamiki
N=120;  % horyzont predykcji
Nu=20;  % horyzont sterowania
U1Max = 80;
U2Max = 100;
dUmax = 0.3;
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
%transmitancja z opoznieniami
G1=tf(sys1);

Gz=c2d(G1,Tp,metoda);
% Dyskretna odpowiedŸ skokowa wielowymiarowa wg konwencjiMatlaba:
Ystep=step(Gz,Tp*D); % macierz o wymiarze (timefinal/Tp+1) x ny x nu na odcinku
Ystep1(:,1:2,1) = Ystep(:,1:2,1);
Ystep1(:,1:2,2) = Ystep(:,1:2,3);
% czasu od t=0 do t=timefinal z krokiem Tp
% Dyskretna odpowiedŸ skokowa macierzowa S o wymiarze ny x nux timefinal/Tp, pierwsza
% macierz dla czasu dyskr.k=1, ostatnia dla czasu dyskr. k=D=timefinal/Tp (D macierzy):
[nt,ny,nu]=size(Ystep1);  % nt=D+1 (Y zawiera te¿ wartoœæ wyjœæ wchwili 0, macierz S nie)
S=zeros(ny,nu,nt-1);
for i=1:ny
for j=1:nu
S(i,j,:)=Ystep1(2:nt,i,j);
end
end
% Wartosci poczatkowe, inicjalizacja macierzy
Yzad = zeros(N*ny,1); 
    Yzad(1,1)=hsp(j);
    Yzad(2,1)=Tvsp(j);
%%% Ustawienie wartoœci zadanej
Y = zeros(N*ny,1); %% wartosc do regulatora
dUp = zeros(2*(D-1),1);
[M,Mp] = DMCmatrices(S,N,Nu);
% Podstawienie do Lambda1
lambda = zeros(Nu*nu,Nu*nu);
wspLambda = 130;
for i=1:nu*Nu
    lambda(i,i) = wspLambda;
end
%%% Wyznaczenie macierzy H wartoœci 
Mpom = (M'*M)+lambda;
H =  2*Mpom;
%%% Wyznaczenie macierzy I 
I = zeros(nu,nu);
for i=1:nu
    I(i,i)=1;
end
%%% Wyznaczenie macierzy J wartoœci
J = zeros(Nu,Nu);
for i=1:2:Nu*nu
    for j=1:2:i+1
    J(i,j) = I(1,1);
    J(i,j+1) = I(1,2);
    J(i+1,j) = I(2,1);
    J(i+1,j+1) = I(2,2);
    end
end
AM = [-J;J;-M;M;];

U = zeros(Tsym/Tp,2);
U(1,1) = Fh0;
U(1,2)= Fc0;

    for i=1:2:N*ny
    Y(i,1)=h0;
    Y(i+1,1)=T0;
    end
    Ymin = zeros(N*ny,1);
    Ymax = zeros(N*ny,1);
      for i=1:2:N*ny
    Ymin(i,1)=0;
    Ymin(i+1,1)=10;
      end
    for i=1:2:N*ny
    Ymax(i,1)=40;
    Ymax(i+1,1)=60;
    end
Upom = zeros(nu*Nu,1);
Umax = zeros(nu*Nu,1);
Umin = zeros(nu*Nu,1);
for i=1:2:nu*Nu
    Umax(i,1)=U1Max;
    Umax(i+1,1)=U2Max;
    Umin(i,1)=0;
    Umin(i+1,1)=0;
end
    
for j=2:Tsym/Tp
     
    for i=1:2:N*ny
    Yzad(i,1)=hsp(j);
    Yzad(i+1,1)=Tvsp(j);
    end
    
    Y0 = Y + Mp * dUp;           
    f = -2*M'*(Yzad-Y0);
    for i=1:2:nu*Nu
    Upom(i,1)=U(j-1,1);
    Upom(i+1,1)=U(j-1,2);
    end
    b = [-Umin+Upom;Umax-Upom;-Ymin+Y0;Ymax-Y0;];

    dU = quadprog(H,f,AM,b,[],[],-dUmax,dUmax);
    U(j,1) = U(j-1,1) + dU(1,1);
    U(j,2) = U(j-1,2) + dU(2,1);
FH(j) = U(j,1);
if(j>tauh/Tp)
Fh_v (j) = FH(j-tauh/Tp);
end
Fc_v (j) = U(j,2);
           h_v(j)=Az(1,1)*h_v(j-1)+Bz(1,1)*Fh_v(j-1)+Bz(1,3)*Fc_v(j-1)+Bz(1,5)*Fd_v(j-1)+Az(1,3)*(-alfa*sqrt(h0)/(A*2));
    h_out(j)=Cz(1,1)*h_v(j)+Dz(1,1)*Fh_v(j)+Dz(1,3)*Fc_v(j)+Dz(1,5)*Fd_v(j)+Cz(1,3)*(-alfa*sqrt(h0)/(A*2));
    T_v(j)=Az(2,1)*h_v(j-1)+Az(2,2)*T_v(j-1)+Bz(2,1)*Fh_v(j-1)+Bz(2,2)*Th_v(j-1)+Bz(2,3)*Fc_v(j-1)+Bz(2,4)*Tc_v(j-1)+Bz(2,5)*Fd_v(j-1)+Bz(2,6)*Td_v(j-1) + Az(2,3)*(-alfa*sqrt(h0)/(A*2));
    if(j>tau/Tp)
        Tout_v(j)=Cz(2,1)*h_v(j-tau/Tp)+Cz(2,2)*T_v(j-tau/Tp)+Dz(2,1)*Fh_v(j-tau/Tp)+Dz(2,2)*Th_v(j-tau/Tp)+Dz(2,3)*Fc_v(j-tau/Tp)+Dz(2,4)*Tc_v(j-tau/Tp)+Dz(2,5)*Fd_v(j-tau/Tp)+Dz(2,6)*Td_v(j-tau/Tp)+Cz(2,3)*(-alfa*sqrt(h0)/(A*2));
    end
    
    dUp(3:2:end,1) = dUp(1:2:end-2,1);
    dUp(4:2:end,1) = dUp(2:2:end-2,1);
    dUp(1,1) = dU(1,1);
    dUp(2,1) = dU(2,1);
    for i=1:2:N*ny
         Y(i,1)=h_out(j);
         Y(i+1,1)=Tout_v(j);
    end
end

figure; plot((0:(length(h_out)-1))*Tp,h_out,'Color','black');
hold on;    plot((0:(length(hsp)-1))*Tp,hsp,'Color','b');
legend('h_o_u_t','hsp','Location','northeast')
title('Wysokosc slupa cieczy')
xlabel('t[s]');
ylabel('h_o_u_t [cm]');

figure; plot((0:(length(Tout_v)-1))*Tp,Tout_v,'Color','black');
hold on;    plot((0:(length(Tvsp)-1))*Tp,Tvsp,'Color','b');
legend('T_o_u_t','Tsp','Location','northeast')
title('Temperatura wyjsciowa')
xlabel('t[s]');
ylabel('T_o_u_t [C]');

figure;  plot((0:(length(Fh_v)-1))*Tp,Fh_v,'Color','blue');
hold on; plot((0:(length(Fc_v)-1))*Tp,Fc_v,'Color','black');
legend('F_h','F_c','Location','northeast')
title('Sterowania')
xlabel('t[s]');
ylabel('F [cm^3/s]');