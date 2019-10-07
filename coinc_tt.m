function [re3ba5tt, re5ba3tt, re3re5tt, ba3ba5tt ] = coinc_tt(pulse_re3, pulse_re5, pba3_vec, pba5_vec)
%This function calculates the coincidence probability in t1-t2 

%Input Variables
%   pulse and background vectors from channel 3 and channel 5

%Output Variables
%   matrix: P(t1,t2) coincidence matrix

n=length(pulse_re3);
re3ba5tt=zeros(n,n);
re5ba3tt=zeros(n,n);
re3re5tt=zeros(n,n);
ba3ba5tt=zeros(n,n);


for i=1:n
    
    for j=1:n
        
        re3ba5tt(i,j)=pulse_re3(i)*pba5_vec(j);
        re5ba3tt(i,j)=pulse_re5(j)*pba3_vec(j);
        re3re5tt(i,j)=pulse_re3(i)*pulse_re5(j);
        ba3ba5tt(i,j)=pba3_vec(i)*pba5_vec(j);
    end
end

% re3ba5tt=pulse_re3.*pba5_vec';
% re5ba3tt=pulse_re5.*pba3_vec';
% re3re5tt=pulse_re3.*pulse_re5';
% ba3ba5tt=pba3_vec.*pba5_vec';

end

