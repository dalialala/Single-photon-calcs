function[vis,  vis_err, par, perp, tau_gate, gate_par, mpar, mperp] = hom_vis(tau_2,tau_3, pt, numer_g3_par, numer_g2_par, numer_g3_perp, numer_g2_perp, winc, file_list_par,file_list_perp, background_flag)
%This function calculates HOM visibility and if background flag=1 it
%calculates the background coincidence per bin and subtracts it
%winc is how wide you want to sum over for coincidences around 0
%win is how wide you want to sum over for background coincidences
%background_flag is a boolean dictating if you want to subtract background

t=tau_2*1e6;
t1=tau_3*1e6;

%Check max tau
max_tau = t(1,end);

%Find t=0
t0=find(t >=0, 1);

%How many cycles there are for t>0
max_bg=floor(max_tau/pt)-1;

%bin time
dt=abs(t(1)-t(2)); 
    
%max_bg=1;
%Find indices for a window (winc) around 0
t0l=find(t >=(-winc), 1)-1;
t0h=find(t >=(winc), 1);
npc=t0h-t0l+1;

t1gate=winc;

HOM=1;
%Define coincidence vector
if dt< 0.05
    tp=4.5785-t1gate;
else
    tp=4.5538-t1gate;
end
[corr_par, single, ~, par] = counts(tau_2, tau_3, numer_g3_par, tp, t1gate, pt ,HOM);
[corr_perp, ~, ~, perp] = counts(tau_2, tau_3, numer_g3_perp, tp, t1gate, pt ,HOM);

% figure
% semilogy(t1,single)

%initialize vectors to find pulse area for |t|>0
index=50;
par_vec=zeros(max_bg-index+1, 2);
perp_vec=zeros(max_bg-index+1, 2);

for i=index:max_bg
        %find relevant time indices for each pulse window
        tpph=find(t >=((i-1)*pt+winc), 1);
        tppl=find(t >=((i-1)*pt-winc), 1);
        tpnl=find(t >=(-(i-1)*pt-winc), 1)-1;
        tpnh=find(t >=(-(i-1)*pt+winc), 1)-1;
        %fill vectors
        par_vec(i-index+1, 1)=sum(par(tppl:tpph));
        par_vec(i-index+1, 2)=sum(par(tpnl:tpnh));
        perp_vec(i-index+1, 1)=sum(perp(tppl:tpph));
        perp_vec(i-index+1, 2)=sum(perp(tpnl:tpnh));        
end

n=2*length(perp_vec);

%normalized area for parallel pulses
%normalized area for perpendicular pulses


area_par=mean(mean(par_vec));
area_perp=mean(mean(perp_vec));


area_par0=sum(par(t0l:t0h));

area_perp0=sum(perp(t0l:t0h));


if background_flag == true
    
    t1gate=winc;
    if dt< 0.05
        tp=4.3978-t1gate;
        tp=4.5785-t1gate;
    else
        tp=4.0304-t1gate;
        tp=4.5538-t1gate;
    end
%     sgate=1;
%     [tau_gate, gate_par, counter_par, pbce_par] = back_profile(tau_2, numer_g2_par, file_list_par, pt, tp, t1gate, sgate, HOM);
%     [tau_gate, gate_perp, counter_perp, pbce_perp] = back_profile(tau_2, numer_g2_perp, file_list_perp, pt, tp, t1gate, sgate, HOM);
% 
%     i1=find(tau_gate >=-winc, 1)-1;
%     i2=find(tau_gate >=winc, 1);
%     
%     bparu=sum(gate_par(i1:i2));
%     bpar_err=pbce_par*(i2-i1+1);
%     
%     bperpu=sum(gate_perp(i1:i2));
%     bperp_err=pbce_perp*(i2-i1+1);
    
    [tau_gate, gate_par, bparu, bpar_err] = int_sc(t, numer_g2_par, pt, tp, winc, file_list_par, HOM);
    [tau_gate, gate_perp, bperpu, bperp_err] = int_sc(t, numer_g2_perp, pt, tp, winc, file_list_perp, HOM);
    
    x=round(200/dt);
    mpar=max(par(t0+x:end));
    %mpar2=max(gate_par);
    mperp=max(perp(t0+x:end));
    
    figure
    semilogy(tau_gate, gate_par/mpar)
    hold
    semilogy(t, perp/mperp)
    semilogy(t, par/mpar)
    xlim([-5, 5])
    ylim([5e-4, 1])
      
    
else
    bparu=0;
    bpar_err=0;
    bperpu=0;
    bperp_err=0;
    tau_gate=0;
    gate_par=0;
    mpar=1;
    mperp=1;
end


%Calculate visibility
num1=area_par0-bparu;
den1=area_par-bparu;
num2=area_perp0-bperpu;
den2=area_perp-bperpu;
ratio1=num1/den1;
ratio2=num2/den2;

vis=1-ratio1/ratio2;

%Calculate error in visibility
n1=sqrt(area_par0+bpar_err^2);
%d1=sqrt(area_par+bpar_err^2);
d1=sqrt(area_par/n+bpar_err^2);
n2=sqrt(area_perp0+bperp_err^2);
%d2=sqrt(area_perp+bperp_err^2);
d2=sqrt(area_perp/n+bperp_err^2);
r1=ratio1*sqrt((n1/num1)^2+(d1/den1)^2);
r2=ratio2*sqrt((n2/num2)^2+(d2/den2)^2);
vis_err=(ratio1/ratio2)*sqrt((r1/ratio1)^2+(r2/ratio2)^2);

end