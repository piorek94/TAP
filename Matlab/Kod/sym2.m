function dyl = sym2( t,y,Fh,Th,Fc,Tc,Fd,Td,alfa,A,Tc0,Th0,Td0,Fc0,Fh0,Fd0,h0,T0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dyl = zeros(2,1);    % a column vector
dyl(1) = (Fh+Fc+Fd-0.5*alfa*sqrt(h0) - 0.5*alfa*y(1)/sqrt(h0))/A;
zm1=(Th0-T0)*Fh+Fh0*Th+(Tc0-T0)*Fc+Fc0*Tc+(Td0-T0)*Fd+Fd0*Td;
dyl(2) = ( ( -( (Fh0*Th0+Fc0*Tc0+Fd0*Td0- (Fh0+Fc0+Fd0) *T0) *y(1) ) /h0) -(Fh0+Fc0+Fd0)*y(2)+zm1)/(A*h0);
end

