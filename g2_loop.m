binv=[105];
vec=[];
err=[];
wind=[];
ts=[];
for i=1:1
tp=(tau_3(174)*1e6-binv(i)*dt);
winc=tau_3(174)*1e6-tp;
t1gate=winc;
[corr, single, t1gate_vec, t1tgate_vec] = counts(tau_2, tau_3, numer_g3, tp, t1gate, pt);
[g20,  g20_err, area_par0, area_par] = g2_pulse(tau, pt, tp, t1tgate_vec, numer_g2, winc, file_list, background_flag);
ts=[ts, tp];
wind=[wind, winc];
vec=[vec, g20];
err=[err, g20_err];
end