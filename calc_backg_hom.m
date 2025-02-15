function [pbc, pbce, pre3, pre3e, pre5, pre5e] = calc_backg_hom(pulse3, pulse5, dt,  tb, tp, counter, nexp)
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

%   Input Variables
%   pbc: Probability of getting a background coincidence 
%   pbce: Uncertainty of probability of getting a background coincidence 
%   pre3: Probability of getting a photon click in channel 3
%   pre3e: Uncertainty of probability of getting a photon click in channel 3
%   pre5: Probability of getting a photon click in channel 3
%   pre5e: Uncertainty of probability of getting a photon click in channel 3


%Probability of singles
pulse3=pulse3/nexp;
pulse5=pulse5/nexp;

%%%%Probability of background clicks and uncertainty
%Background window start index
i1=round(tb/dt);
%Background window end index
i2=round((tb+1.2)/dt)

%Bacground sum in window
pba3=sum(pulse3(i1:i2));
pba5=sum(pulse5(i1:i2));
%Background window bin number
nbw=abs(i2-i1+1);
%Actual background probability and uncertainty in a bin dt
pba3e=sqrt(pba3)/sqrt(nexp*counter)/nbw;
pba5e=sqrt(pba5)/sqrt(nexp*counter)/nbw;
pba3=pba3/nbw
pba5=pba5/nbw


pulse_re3=pulse3-pba3;
pulse_re5=pulse5-pba5;
%%%%Probability of background clicks and uncertainty
%Pulse window start index
i1=round(tp/dt);
%Pulse window end index
i2=round((tp+2.25)/dt);
%Actual real  probability and uncertainty in a bin dt
pre3=sum(pulse_re3(i1:i2));
pre5=sum(pulse_re5(i1:i2));
pre3e=sqrt(pre3/nexp/counter);
pre5e=sqrt(pre5/nexp/counter);

%Probability of background coincidence
pbc=(pre3*pba5+pre5*pba3+pba3*pba5);
%Probability of background coincidence uncertainty
un1=pre3*pba5*sqrt((pre3e/pre3)^2+(pba5e/pba5)^2);
un2=pre5*pba3*sqrt((pre5e/pre5)^2+(pba3e/pba3)^2);
un3=pba5*pba3*sqrt((pba5e/pba5)^2+(pba3e/pba3)^2);
pbce=sqrt(un1^2+un2^2+un3^2);

end
