function [back_3, back_5, pba3, pba3e, pba5,  pba5e, pbce, pulse_re3, pulse_re5, pre3, pre3e, pre5, pre5e] = calc_backg(pulse3, pulse5, dt,  tb, tp, t1gate, counter, nexp, HOM)
%This fucntion calculates the "real" probability of getting a photon click
%and the probability of getting a background coincidence in a bin "dt" from a vector which
%contains the single counts over a cycle of duration pt and averaged
%"n=counter" times.

%   Input Variables
%   pulse3: Averaged single counts from channel 3 
%   pulse5: Averaged single counts from channel 5 
%   dt: bin time
%   tb: Background window starts
%   tp: Pulse window starts
%   counter: numer of times the cycle has been averaged
%   tot_files: total numer of experiments
%   HOM: Boolean

%   Input Variables
%   pba3: Probability of getting a background click in channel 3
%   pba3e: Uncertainty of probability of getting a background click in channel 3
%   pba5: Probability of getting a background click in channel 3
%   pba5e: Uncertainty of probability of getting a background click in channel 3
%   pbce: Unceratinty of background coincedence in a bin dt 
%   pulse_re3: single counts from channel 3 with background subtraction
%   pulse_re5: single counts from channel 5 with background subtraction
%   pre3: Probability of getting a photon click in channel 3
%   pre3e: Uncertainty of probability of getting a photon click in channel 3
%   pre5: Probability of getting a photon click in channel 3
%   pre5e: Uncertainty of probability of getting a photon click in channel 3


%Probability of singles
pulse3=pulse3/nexp;
pulse5=pulse5/nexp;

%%%%Probability of background clicks and uncertainty
%Background window start index
if HOM == true
    if dt< 0.05
        ib1=18;
    else
        ib1=1;
    end
else
    ib1=round(tb/dt);
    ib1=210;
    %ib1=35
end

%ib1=210;
%Background window end index
if HOM == true
    if dt<0.05
        ib2=88;
    else
        ib2=34;
    end
else
    ib2=40;
    ib2=length(pulse3);
    %ib2=100
    
end

%Bacground sum in window
pba3=sum(pulse3(ib1:ib2));
pba5=sum(pulse5(ib1:ib2));
%Background window bin number
nbw=ib2-ib1+1;
%Actual background probability and uncertainty in a bin dt
pba3e=sqrt(pba3)/sqrt(nexp*counter)/nbw;
pba5e=sqrt(pba5)/sqrt(nexp*counter)/nbw;
pba3=pba3/nbw;
pba5=pba5/nbw;


pulse_re3=pulse3-pba3;
pulse_re5=pulse5-pba5;

%%%%Probability of photon clicks and uncertainty
%Pulse window start index
i1=round(tp/dt)+1;
%Pulse window end index
i2=round((tp+t1gate)/dt)+1;
%Actual real  probability and uncertainty in a bin dt
pre3=sum(pulse_re3(i1:i2));
pre5=sum(pulse_re5(i1:i2));
pre3e=sqrt(pre3/nexp/counter);
pre5e=sqrt(pre5/nexp/counter);

%Leakage window start index
tli=3.69;
il1=round(tli/dt);
%Leakage window end index
tle=tli+0.56;
il2=round(tle/dt);

% back_3=pulse3;
% back_5=pulse5;
% back_3(i1:i2)=pba3;
% back_3(il1:il2)=pba3;
% back_5(i1:i2)=pba5;
% back_5(il1:il2)=pba5;
back_3=pba3*ones(1, length(pulse3));
back_5=pba5*ones(1, length(pulse3));



%Probability of background coincidence per bin
pbc=(pre3*pba5+pre5*pba3+pba3*pba5);
%Probability of background coincidence uncertainty
un1=pre3*pba5*sqrt((pre3e/pre3)^2+(pba5e/pba5)^2);
un2=pre5*pba3*sqrt((pre5e/pre5)^2+(pba3e/pba3)^2);
un3=pba5*pba3*sqrt((pba5e/pba5)^2+(pba3e/pba3)^2);
pbce=sqrt(un1^2+un2^2+un3^2);

end
