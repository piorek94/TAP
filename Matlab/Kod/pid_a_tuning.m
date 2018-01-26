
%nastawy regulatorow
kp_h = 0.514;
ki_h = 0.0013;
kp_T = 0.298;
ki_T = 0.0011;

kp_h_lower = 0.1;
kp_h_upper = 3;
ki_h_lower = 0.0005;
ki_h_upper = 0.01;

kp_T_lower = 0.1;
kp_T_upper = 3;
ki_T_lower = 0.0005;
ki_T_upper = 0.01;

N = 4;
testPoints = N^4;
%kp_h ki_h kp_T ki_T
%parameters = (kp_h_lower : (kp_h_upper-kp_h_lower)/(testPoints-1) : kp_h_upper)';
%parameters = [parameters (ki_h_lower : (ki_h_upper-ki_h_lower)/(testPoints-1) : ki_h_upper)'];
%parameters = (kp_T_lower : (kp_T_upper-kp_T_lower)/(testPoints-1) : kp_T_upper)';
%parameters = [parameters (ki_T_lower : (ki_T_upper-ki_T_lower)/(testPoints-1) : ki_T_upper)'];


sweep = 1;
if sweep == 1
    % sweep PIDs parameters
    %kp_h ki_h kp_T ki_T
    parameters = zeros(testPoints, 5);
    i = 1;
    for kp_h = (kp_h_lower : (kp_h_upper-kp_h_lower)/(N-1) : kp_h_upper)
       for ki_h = (ki_h_lower : (ki_h_upper-ki_h_lower)/(N-1) : ki_h_upper)
          for kp_T = (kp_T_lower : (kp_T_upper-kp_T_lower)/(N-1) : kp_T_upper)
              for ki_T = (ki_T_lower : (ki_T_upper-ki_T_lower)/(N-1) : ki_T_upper)
                wsk = BLTwsklog(G1, kp_h, ki_h, kp_T, ki_T, 0);
                parameters(i,:) = [kp_h, ki_h, kp_T, ki_T, wsk];
                disp(i);
                i = i + 1;
              end
          end
       end
    end
end

% find parameters resulting in correct BLT value
bestParameters = zeros(1,5);
for i = (1 : testPoints)
    wsk = parameters(i,5);
    if(abs(wsk - 4.0) < abs(bestParameters(5) - 4.0))
        bestParameters = parameters(i,:)
    end
end

BLTwsklog(G1, bestParameters(1), bestParameters(2), bestParameters(3), bestParameters(4), 1);