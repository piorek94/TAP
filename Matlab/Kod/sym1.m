function dy = sym1( t,y,Fh,Th,Fc,Tc,Fd,Td,alfa,A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dy = zeros(2,1);    % a column vector
dy(1) = (Fh+Fc+Fd-alfa*sqrt(y(1)))/A;
dy(2) = (Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*y(2))/(A*y(1));
end

