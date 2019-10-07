function [tau_gate, gate_corr] = back_profile(tau, numer_g2, file_list, pt, tp, t1gate, sgate)

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
tn=round(tend/pt);

%Flip time in singles
singles_3 = double(fliplr(numer_g2(2,:)));
singles_5 = double(fliplr(numer_g2(3,:)));

singles_3 = double((numer_g2(2,:)));
singles_5 = double((numer_g2(3,:)));


[pulse_3 ,counter] = average_cycle(tn, t0, dt, pt, singles_3);
[pulse_5 ,counter] = average_cycle(tn, t0, dt, pt, singles_5);

tpulse=t(t0:t0+round(pt/dt)-1);
tpulse2=t(t0:t0+round(2*pt/dt)-1);
tpulse3=t(t0:t0+round(3*pt/dt)-1);


% tp=1.46;
% tb=3.75;
% tp=1.41;
tb=dt;

%Calculates background probability on both channels
[pba3vec, pba5vec, pba3, pba3e, pba5, pba5e, pulsere3, pulsere5, pre3, pre3e, pre5, pre5e] = calc_backg(pulse_3, pulse_5, dt,  tb, tp, counter, nexp);

gatet=t1gate;

n=length(pulsere3);
pulse_re3=zeros(1,n);
pulse_re5=zeros(1,n);
pba3_vec=zeros(1,n);
pba5_vec=zeros(1,n);

if sgate == 1
    itp=ceil(tp/dt)+1;
    itg=ceil((tp+gatet)/dt)+1;
    pulse_re3(itp: itg)=pulsere3(itp: itg);
    pulse_re5(itp: itg)=pulsere5(itp: itg);
    pba3_vec(itp: itg)=pba3vec(itp: itg);
    pba5_vec(itp: itg)=pba5vec(itp: itg);
    tau1=t(t0:t0+round(pt/dt));
    round(pt/dt)+1;
    tau1(itp)
    tau1(itg)
else
    pulse_re3=pulsere3;
    pulse_re5=pulsere5;
    pba3_vec=pba3vec;
    pba5_vec=pba5vec;
    gatet=2.0;  
    itp=1;
    itg=round(pt/dt);
end


%Background coincidences in t1 t2 space
 n=length(pulse_re3);

[re3ba5tt, re5ba3tt, re3re5tt, ba3ba5tt ] = coinc_tt(pulse_re3, pulse_re5, pba3_vec, pba5_vec);

[re3ba5] = coinc_tau(tpulse3, [re3ba5tt, re3ba5tt, re3ba5tt]);
[re5ba3] = coinc_tau(tpulse3, [re5ba3tt, re5ba3tt, re5ba3tt]);
[re3re5] = coinc_tau(tpulse3, [re3re5tt, re3re5tt, re3re5tt] );
[ba3ba5] = coinc_tau(tpulse3, [ba3ba5tt, ba3ba5tt, ba3ba5tt]);

tott=((re3ba5+re5ba3+re3re5+ba3ba5)*nexp);

figure
[T1, Tau]=meshgrid(tpulse3, tpulse);
surf(T1, Tau, tott), colorbar, view(2), EdgeColor = 'none';
shading interp

tau_vec=t(t0-round(1.5*pt/dt):t0+round(1.5*pt/dt));

figure
semilogy(sum(tott(1:end,:), 1))


pulse_re30=pulse_re3;
pulse_re50=pulse_re5;

itp1=ceil(tp/dt)+1;
itg1=ceil((tp+2.0)/dt);

pulse_re30(itp1:itg1)=0;
pulse_re50(itp1:itg1)=0;

pulse_re5t=[pulse_re5, pulse_re50, pulse_re5];
pulse_re3t=[pulse_re3, pulse_re30, pulse_re3];

pba3_vec=[pba3_vec, pba3_vec, pba3_vec];
pba5_vec=[pba5_vec, pba5_vec, pba5_vec];


[re3ba5t0, re5ba3t0, re3re5t0, ba3ba5t0] = coinc_tt(pulse_re3t, pulse_re5t, pba3_vec, pba5_vec);

[re3ba5_0] = coinc_tau(tpulse3, re3ba5t0);
[re5ba3_0] = coinc_tau(tpulse3, re5ba3t0);
[re3re5_0] = coinc_tau(tpulse3, re3re5t0);
[ba3ba5_0] = coinc_tau(tpulse3, ba3ba5t0);

tot=((re3ba5_0+re5ba3_0+re3re5_0+ba3ba5_0)*nexp);
%tot0=round((re3ba5+re5ba3+ba3ba5)*nexp);

tot1=tot(2*n+1:end, :);
tot2=tot(n+1:2*n, :);

ix1=round(0.72/dt);
ix2=round(1/dt);
iy1=round(7/dt);
iy2=round(6.1/dt);
tot2(ix1:end, iy1:end)=0;
tot2(1:ix2, iy2:end)=0;
tot=tot1+ tot2;

tot=[tot(:, round(2.5*n):end), tot(:, 1:round(2.5*n)-1)];

[T1, Tau]=meshgrid(tpulse3-5, tpulse);
figure
surf(T1, Tau,tot), colorbar, view(2), EdgeColor = 'none';
shading interp

gate_corr=ceil(sum(tot(itp:itg,:), 1));

figure
tau_gate=linspace(-1.5*pt, 1.5*pt, 3*n);
semilogy(tau_gate, gate_corr)
end

