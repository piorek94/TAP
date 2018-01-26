%clear
clc
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
tau=0;
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
dUmax = 1;
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
Yzad = zeros(N,2); 
%%% Ustawienie wartoœci zadanej
for i=1:N
    Yzad(i,1)=18.65;
    Yzad(i,2)=33.16;
end
Y = zeros(N,2); %% wartosc do regulatora
dUp = zeros(D-1,2);
U = zeros(D,2);
[M,Mp] = DMCmatrices(S,N,Nu);
% Podstawienie do Lambda1
lambda = zeros(Nu*nu,Nu*nu);
wspLambda1 = 130;
wspLambda2 = 130;
for i=1:2:nu*Nu
    lambda(i,i) = wspLambda1;
    lambda(i+1,i+1) = wspLambda2;
end
WspFi = zeros(N*ny,N*ny);
%%% Wyznaczenie macierzy K wartoœci
for i=1:2:ny*N
    WspFi(i,i) = 1;
    WspFi(i+1,i+1) = 1;
end
Mpom = (M'*WspFi*M)+lambda;
K = Mpom\M'; 

K1 = zeros(1,N*2);
for i=1:N*2
K1(1,i) = K(1,i);
K1(2,i) = K(2,i);
end
Ke = zeros(2,2);
for i=1:N
Ke(1,1) = Ke(1,1)+K1(1,i*2-1);
Ke(1,2) = Ke(1,2)+K1(1,i*2);
Ke(2,1) = Ke(2,1)+K1(2,i*2-1);
Ke(2,2) = Ke(2,2)+K1(2,i*2);
end
Mpj = zeros(2*N,2,D-1);
for i=1:D-1
Mpj(:,:,i) = Mp(:,2*i-1:2*i);
end
Ku = zeros(ny,nu,D-1);
for i=1:D-1
Ku(:,:,i) = K1*Mpj(:,:,i);
end
%%%%%%%%%%%%%%%%%%%TUTAJ SKONCZYLEM - ogarnac Ku - odpowiednio odczytac Mpj
u = zeros(D,2);
for i=1:D
u(i,1) = Fh0;
u(i,2)= Fc0;
end
y = zeros(D,2);
for i=1:D
y(i,1) = h0;
y(i,2)= T0;
end
e = zeros(2,1);
for j=2:Tsym/Tp
     
        Yzad(1,1) = hsp(j);
        Yzad(1,2) = Tvsp(j);
          e(1,1) = Yzad(1,1) - y(j-1,1);
          e(2,1) = Yzad(1,2) - y(j-1,2);
            
   upom = u;
        for i=1:D-1
         u(D-i,1)=upom(D-i+1,1);
         u(D-i,2)=upom(D-i+1,2);
        end
        for i=1:D-1
        Kupom(:,i) = Ku(:,:,i)*dUp(i,:)';
        end
        EKu = zeros(2,1);
        for i=1:D-1
            EKu(1,1) = EKu(1,1)+Kupom(1,i);
            EKu(2,1) = EKu(2,1)+Kupom(2,i);
        end
        Kepom = Ke*e;
        %%% Obliczenie trajektorii swobodnej Y0(k)
        dU(1,1) = Kepom(1,1) - EKu(1,1);
        dU(2,1) = Kepom(2,1) - EKu(2,1);
        %%%Ograniczenia na przyrost sterowania
        if(dU(1,1)>dUmax)
            dU(1,1) = dUmax;
        end
        if(dU(1,1)<-dUmax)
            dU(1,1) = -dUmax;
        end
        if(dU(2,1)>dUmax)
            dU(2,1) = dUmax;
        end
        if(dU(2,1)<-dUmax)
            dU(2,1) = -dUmax;
        end
        %%% Wyliczenie wyjscia - sygnalu sterowania
u(D,1) = u(D-1,1) + Kepom(1,1) - EKu(1,1);
u(D,2) = u(D-1,2) + Kepom(2,1) - EKu(2,1);
%%% Rzutowanie sterowania na ograniczenia
if(u(D,1)<0)
    u(D,1) = 0;
end
if(u(D,1)>U1Max)
    u(D,1) = U1Max;
end
if(u(D,2)<0)
    u(D,2) = 0;
end
if(u(D,2)>U2Max)
    u(D,2) = U2Max;
end
FH(j) = u(D,1);
if(j>tauh/Tp)
Fh_v (j) = FH(j-tauh/Tp);
end
Fc_v (j) = u(D,2);
           h_v(j)=Az(1,1)*h_v(j-1)+Bz(1,1)*Fh_v(j-1)+Bz(1,3)*Fc_v(j-1)+Bz(1,5)*Fd_v(j-1)+Az(1,3)*(-alfa*sqrt(h0)/(A*2));
    h_out(j)=Cz(1,1)*h_v(j)+Dz(1,1)*Fh_v(j)+Dz(1,3)*Fc_v(j)+Dz(1,5)*Fd_v(j)+Cz(1,3)*(-alfa*sqrt(h0)/(A*2));
    y(j,1) = h_out(j);
    T_v(j)=Az(2,1)*h_v(j-1)+Az(2,2)*T_v(j-1)+Bz(2,1)*Fh_v(j-1)+Bz(2,2)*Th_v(j-1)+Bz(2,3)*Fc_v(j-1)+Bz(2,4)*Tc_v(j-1)+Bz(2,5)*Fd_v(j-1)+Bz(2,6)*Td_v(j-1) + Az(2,3)*(-alfa*sqrt(h0)/(A*2));
    if(j>tau/Tp)
        Tout_v(j)=Cz(2,1)*h_v(j-tau/Tp)+Cz(2,2)*T_v(j-tau/Tp)+Dz(2,1)*Fh_v(j-tau/Tp)+Dz(2,2)*Th_v(j-tau/Tp)+Dz(2,3)*Fc_v(j-tau/Tp)+Dz(2,4)*Tc_v(j-tau/Tp)+Dz(2,5)*Fd_v(j-tau/Tp)+Dz(2,6)*Td_v(j-tau/Tp)+Cz(2,3)*(-alfa*sqrt(h0)/(A*2));
        y(j,2) = Tout_v(j);
    end
    dUp(2:end,1) = dUp(1:end-1,1);
    dUp(2:end,2) = dUp(1:end-1,2);
    dUp(1,1) = dU(1,1);
    dUp(1,2) = dU(2,1);
end

figure; plot((0:(length(h_out)-1))*Tp,h_out,'Color','black','LineWidth',1.5);
hold on;    plot((0:(length(hsp)-1))*Tp,hsp,'Color','b','LineWidth',1.5);
hold on;    plot((0:(length(h_out)-1))*Tp,h_outOp,'Color','red','LineWidth',1.5);
% hold on;    plot((0:(length(h_out)-1))*Tp,h_outN60,'Color','green','LineWidth',1.5);
% hold on;    plot((0:(length(h_out)-1))*Tp,h_outNu3,'Color','magenta','LineWidth',1.5);
legend('Obiekt bez opóŸnienia tau = 0','h_s_p','Obiekt z opóŸnieniem tau = 270s','Location','northeast')
title('Wysokosc slupa cieczy','FontSize',12);
xlabel('t[s]','FontSize',12);
ylabel('h_o_u_t [cm]','FontSize',12);

figure; plot((0:(length(Tout_v)-1))*Tp,Tout_v,'Color','black','LineWidth',1.5);
hold on;    plot((0:(length(Tvsp)-1))*Tp,Tvsp,'Color','b','LineWidth',1.5);
hold on;    plot((0:(length(Tout_v)-1))*Tp,Tout_vOp,'Color','red','LineWidth',1.5);
%hold on;    plot((0:(length(Tout_v)-1))*Tp,Tout_vN60,'Color','green','LineWidth',1.5);
%hold on;    plot((0:(length(Tout_v)-1))*Tp,Tout_vNu3,'Color','magenta','LineWidth',1.5);
legend('Obiekt bez opóŸnienia tau = 0','T_s_p','Obiekt z opóŸnieniem tau = 270s','Location','northeast')
title('Temperatura wyjsciowa','FontSize',12)
xlabel('t[s]','FontSize',12);
ylabel('T_o_u_t [C]','FontSize',12);

figure;  plot((0:(length(Fh_v)-1))*Tp,Fh_v,'Color','blue','LineWidth',1.5);
hold on; plot((0:(length(Fc_v)-1))*Tp,Fc_v,'Color','black','LineWidth',1.5);
hold on; plot((0:(length(Fh_v)-1))*Tp,Fh_vOp,'Color','green','LineWidth',1.5);
hold on; plot((0:(length(Fc_v)-1))*Tp,Fc_vOp,'Color','red','LineWidth',1.5);
legend('F_h bez opóŸnienia tau','F_c bez opóŸnienia tau','F_h z opóŸnieniem tau','F_c z opóŸnieniem tau','Location','northeast')
title('Sterowania','FontSize',12)
xlabel('t[s]','FontSize',12);
ylabel('F [cm^3/s]','FontSize',12);

%figure; plot((0:(length(Ystep1(:,1,1))-1))*Tp,Ystep1(:,1,1),'Color','black','LineWidth',1.5);
%hold on;    plot((0:(length(Ystep1(:,2,1))-1))*Tp,Ystep1(:,2,1),'Color','b','LineWidth',1.5);
%legend('h_o_u_t','T_o_u_t','Location','northeast')
%title('OdpowiedŸ skokowa na skok F_H_i_n','FontSize',12);
%xlabel('t[s]','FontSize',12);
%ylabel('h_o_u_t , T_o_u_t','FontSize',12);

%figure; plot((0:(length(Ystep1(:,1,2))-1))*Tp,Ystep1(:,1,2),'Color','black','LineWidth',1.5);
%hold on;    plot((0:(length(Ystep1(:,2,2))-1))*Tp,Ystep1(:,2,2),'Color','b','LineWidth',1.5);
%legend('h_o_u_t','T_o_u_t','Location','northeast')
%title('Odpowiedz skokowa na skok F_C','FontSize',12)
%xlabel('t[s]','FontSize',12);
%ylabel('h_o_u_t , T_o_u_t','FontSize',12);