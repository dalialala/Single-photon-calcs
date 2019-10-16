binv=[60,  70];
vec=[];
err=[];
wind=[];
ts=[];
for i=2:2
dt=(tau(2)-tau(1))*1e6;
pt=round(5/dt)*dt;
tp=(tau_3(174)*1e6-binv(i)*dt);
winc=binv(i)*dt;
t1gate=winc;
background_flag=1;
[corr, single, t1gate_vec, t1tgate_vec] = counts(tau_2, tau_3, numer_g3, tp, t1gate, pt,0);
[g20,  g20_err, area_par0, area_par] = g2_pulse(tau, pt, tp, t1tgate_vec, numer_g2, winc, file_list, background_flag);
ts=[ts, tp];
wind=[wind, winc];
vec=[vec, g20];
err=[err, g20_err];
end