function[g20,  g20_err, area_par0, area_par] = g2_pulse(tau, pt, tp, vec, numer_g2, winc, file_list, background_flag)
%This function calculates HOM visibility and if background flag=1 it
%calculates the background coincidence per bin and subtracts it
%winc is how wide you want to sum over for coincidences around 0
%win is how wide you want to sum over for background coincidences
%background_flag is a boolean dictating if you want to subtract background

t=tau*1e6;

%Check max tau
max_tau = t(1,end);

%Find t=0
t0=find(t >=0, 1);

%How many cycles there are for t>0
max_bg=floor(max_tau/pt);

%bin time
dt=abs(t(1)-t(2)); 
    
%max_bg=1;
%Find indices for a window (winc) around 0
t0l=find(t >=(-winc), 1)-1;
t0h=find(t >=(winc), 1);
npc=t0h-t0l+1;

%Define coincidence vector
par=double(vec);

%initialize vectors to find pulse area for |t|>0
index=50;
p_vec=zeros(max_bg-index+1, 2);

for i=index:max_bg
        %find relevant time indices for each pulse window
        tpph=find(t >=((i-1)*pt+winc), 1);
        tppl=find(t >=((i-1)*pt-winc), 1);
        tpnl=find(t >=(-(i-1)*pt-winc), 1)-1;
        tpnh=find(t >=(-(i-1)*pt+winc), 1)-1;
        %fill vectors
        p_vec(i-index+1, 1)=sum(par(tppl:tpph));
        p_vec(i-index+1, 2)=sum(par(tpnl:tpnh));
        
end

%normalized area for parallel pulses
area_par=mean(mean(p_vec));
n=2*length(p_vec);

area_par0=sum(par(t0l:t0h));

if background_flag == true
    
    %tp=1.46;
    sgate=1;
    t1gate=winc;
    HOM=0;
%     [tau_gate, gate_corr, counter, pbce] = back_profile(tau, numer_g2, file_list, pt, tp, t1gate, sgate, HOM);
%     i1=find(tau_gate <=-winc, 1, 'last');
%     i2=find(tau_gate <=winc, 1, 'last');
%     
%     bparu=sum(gate_corr(i1:i2));
%     bpar_err=pbce*(i2-i1+1)

    [tau_gate, gate_corr, bparu, bpar_err] = int_sc(t, numer_g2, pt, tp,winc, file_list, HOM);

    
    figure
    semilogy(t, vec)
    hold
    semilogy(tau_gate, (gate_corr))
    xlim([-7, 7])
    ylim([20, 50000])
    
else
    bparu=0;
    bpar_err=0;
    tau_gate=t;
    gate_corr=vec(t0l:t0h);
end


%Calculate g20

num1=area_par0-bparu;
den1=area_par-bparu;

ratio1=num1/den1;

g20=ratio1;

%Calculate error in g20
n1=sqrt(area_par0+bpar_err^2);
d1=sqrt(area_par/n+bpar_err^2);
r1=ratio1*sqrt((n1/num1)^2+(d1/den1)^2);
g20_err=abs(r1);

end