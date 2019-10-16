function [tau_gate, gate_corr, counter, pbce] = back_profile(tau, numer_g2, file_list, pt, tp, t1gate, sgate, HOM)

%Convert time vector to us
t=tau*1e6;

%Bin size
dt=t(2)-t(1);

%tau=0 index
t0=find(t >=0, 1);
  
%Maximum time
tend=t(end);

%Total measurements
tot_files=length(file_list);

%Cycles in 0.6ms with pulse time pt
cycles=0.6e6/pt;

%Number of total experiments
nexp=tot_files*cycles;

%Repetition times for positive times
tn=round(tend/pt)-1;

singles_3 = double((numer_g2(2,:)));

singles_5 = double((numer_g2(3,:)));

%For HOM meas
if HOM == true
    t0a=t0+10;
else
    t0a=t0;
end
    
[pulse_3 ,counter] = average_cycle(t, tn, t0a, dt, pt, singles_3);
[pulse_5 ,counter] = average_cycle(t, tn, t0a, dt, pt, singles_5);


t1 = find(t<=pt,1, 'last')-1;
tpulse=t(t0:t1);
length(tpulse);
t3 = find(t<=3*pt, 1, 'last')-1;
tpulse3=t(t0:t3);
% length(tpulse3)
% length(pulse_3)

tb=dt;

%Calculates background probability on both channels
[pba3vec, pba5vec, pba3, pba3e, pba5, pba5e, pbce, pulsere3, pulsere5, pre3, pre3e, pre5, pre5e] = calc_backg(pulse_3, pulse_5, dt,  tb, tp, t1gate, counter, nexp, HOM);

gatet=t1gate;

n=length(pulsere3);
pulse_re3=zeros(1,n);
pulse_re5=zeros(1,n);
pba3_vec=zeros(1,n);
pba5_vec=zeros(1,n);
% % 
% figure
% semilogy(tpulse, pulse_3)

if sgate == 1
    itp=find(tpulse<=tp,1, 'last')+1;
    itg=find(tpulse<=(tp+gatet),1, 'last');
    tpulse(itp)
    tpulse(itg)
    tpulse(itg)-tpulse(itp)
    pulse_re3(itp: itg)=pulsere3(itp: itg);
    pulse_re5(itp: itg)=pulsere5(itp: itg);
    pba3_vec(itp: itg)=pba3vec(itp: itg);
    pba5_vec(itp: itg)=pba5vec(itp: itg);

else
    pulse_re3=pulsere3;
    pulse_re5=pulsere5;
    pba3_vec=pba3vec;
    pba5_vec=pba5vec;
    %gatet=2.0; 
    tp=1.4458;
    gatet=2.2089;
    itp=1;
    itg=n;
end

% figure
% semilogy(tpulse, pulse_re3)
% hold
% semilogy(tpulse, pba3_vec)
% semilogy(tpulse, pulse_re5)
% semilogy(tpulse, pba5_vec)

%Background and real coincidences in t1 t2 space
% [re3ba5tt, re5ba3tt, re3re5tt, ba3ba5tt ] = coinc_tt(pulse_re3, pulse_re5, pba3_vec, pba5_vec);
% 
% [re3ba5] = coinc_tau(tpulse3, [re3ba5tt, re3ba5tt, re3ba5tt]);
% [re5ba3] = coinc_tau(tpulse3, [re5ba3tt, re5ba3tt, re5ba3tt]);
% [re3re5] = coinc_tau(tpulse3, [re3re5tt, re3re5tt, re3re5tt] );
% [ba3ba5] = coinc_tau(tpulse3, [ba3ba5tt, ba3ba5tt, ba3ba5tt]);
% 
% tott=((re3ba5+re5ba3+re3re5+ba3ba5)*nexp);
% 
% figure
% [T1, Tau]=meshgrid(tpulse3, tpulse);
% surf(T1, Tau, tott), colorbar, view(2), EdgeColor = 'none';
% shading interp

% figure
% semilogy(sum(tott(1:end,:), 1))
 
pulse_re30=pulse_re3;
pulse_re50=pulse_re5;
% 
% figure
% semilogy(tpulse, pulse_re3)

itpt=find(tpulse<=tp,1, 'last');
itgt=find(tpulse<=(tp+gatet),1, 'last');

pulse_re30(itpt:itgt)=0;
pulse_re50(itpt:itgt)=0;

pulse_re5t=[pulse_re5, pulse_re50, pulse_re5];
pulse_re3t=[pulse_re3, pulse_re30, pulse_re3];

pba3_vec=[pba3_vec, pba3_vec, pba3_vec];
pba5_vec=[pba5_vec, pba5_vec, pba5_vec];

% figure
% semilogy(tpulse3, pulse_re3t)
% hold
% semilogy(tpulse3, pulse_re5t)
% semilogy(tpulse3, pba3_vec)
% semilogy(tpulse3, pba5_vec)


[re3ba5t0, re5ba3t0, re3re5t0, ba3ba5t0] = coinc_tt(pulse_re3t, pulse_re5t, pba3_vec, pba5_vec);

[re3ba5_0] = coinc_tau(tpulse3, re3ba5t0);
[re5ba3_0] = coinc_tau(tpulse3, re5ba3t0);
[re3re5_0] = coinc_tau(tpulse3, re3re5t0);
[ba3ba5_0] = coinc_tau(tpulse3, ba3ba5t0);



tot=((re3ba5_0+re5ba3_0+re3re5_0+ba3ba5_0)*nexp);
tot1=tot(2*n+1:end, :);
tot2=tot(n+1:2*n, :);
 
ix1=round(0.5/dt);
ix2=round(1/dt);
iy1=round(7/dt);
iy2=round(5/dt);
tot2(ix1:end, iy1:end)=0;
tot2(1:ix2, iy2:end)=0;
tot=tot1+ tot2;

tot=[tot(:, round(2.5*n):end), tot(:, 1:round(2.5*n)-1)];

tot(1:itp-1,:)=0;
tot(itg+1:end, :)=0;
tau_gate=t(t0-floor(1.5*n)-1:t0+floor(1.5*n)-1);


 [T1, Tau]=meshgrid(tau_gate, tpulse);
%  figure
%  surf(T1, Tau,tot1), colorbar, view(2);
%  shading interp
%  
%   figure
%  surf(T1, Tau,tot2), colorbar, view(2);
%  shading interp

gate_corr=(sum(tot(1:end,:), 1));
length(gate_corr);
length(tau_gate);

 figure
 surf(T1, Tau,tot), view(2);
 shading interp
 xlim([-3,3])
 ylim([1.5, 4.0])
 caxis([0 30])
 pbaspect([1 1 1])
 
end

