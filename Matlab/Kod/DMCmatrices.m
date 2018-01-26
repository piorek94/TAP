
function [M,MP]=DMCmatrices(S,N,Nu)
% script "DMCmatrices" macierze M i MP regulatora DMC,
% DANE WEJŒCIOWE (wymagane):
% macierz skoñczonych odpowiedzi skokowych S(ny,nu,D), gdzie:
% ny = dim y; nu = dim u; D - horyzont dynamiki od k=1 do k=D
% N - horyzont predykcji;
% Nu - horyzont sterowania
[ny,nu,D]=size(S); % D - liczba dyskretnych chwil czasu odpowiedzi skokowej
% Macierz dynamiczna M:
M=zeros(N*ny,Nu*nu);
for i=1:N, M((i-1)*ny+1:(i-1)*ny+ny,1:nu)=S(:,:,min(i,D)); end
for i=2:Nu
M(:,(i-1)*nu+1:(i-1)*nu+nu)=[zeros((i-1)*ny,nu); M(1:(N-i+1)*ny,1:nu)];
end
% macierz MP:
MP=zeros(ny*N,nu*(D-1));
for i=1:D-1
for j=1:N
MP((j-1)*ny+1:j*ny,(i-1)*nu+1:i*nu)=S(:,:,min(i+j,D))-S(:,:,i);
end
end