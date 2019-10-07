function[g20,  g20_err, area_par0, area_par] = g2_pulse(tau, pt, numer_g2, winc, file_list, background_flag)
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
par=double(numer_g2(1,:));

%initialize vectors to find pulse area for |t|>0
index=5;
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

npar=1;

%Calculate mean area for pulses with t>pt above certain threshold
th=0.975;
count_p=0;
area_p=0;
for i=1:length(p_vec)
    if p_vec(i,1)>=max(max(p_vec))*th
        area_p=area_p+p_vec(i,1);
        p_vec(i,1);
        count_p=count_p+1; 
    end
    if p_vec(i,2)>=max(max(p_vec))*th
        area_p=area_p+p_vec(i,2);
        p_vec(i,2);
        count_p=count_p+1; 
    end
   
end

%normalized area for parallel pulses
area_par=area_p/count_p;

area_par0=sum(par(t0l:t0h));

if background_flag == true
    
    %Cycles from Setlist
    pt=round(pt, 3);
    cycles=floor(0.6e6/pt);

    
    %total files
    tot_files=length(file_list);
    
    %Nunmber of experiments
    nexp=tot_files*cycles
    
    %Single counts from each channel
    singles3_p=double(fliplr(numer_g2(2,:)));
    singles5_p=double(fliplr(numer_g2(3,:)));
    
    %Get averaged cycle for each channel
    [av_par3 ,counter] = average_cycle(max_bg, t0, dt, pt,  singles3_p);
    [av_par5 ,counter] = average_cycle(max_bg, t0, dt, pt,  singles5_p);
    
    %Background window start time
    tb=3.75;
    
    %Pulse window start time
    tp=1.46;
    
    %Probability of background coincidence and uncertainty
    [pbc, pbce, ~, ~, ~, ~] = calc_backg_hom(av_par3, av_par5, dt,  tb, tp, counter, nexp);
    
    %Background counts in pulse window  
    bparu=npc*pbc*nexp;
    
    %Uncertainty on background coincidence
    bpar_err=npc*pbce*nexp;
    
    
else
    bparu=0;
   
    bpar_err=0;
end


%Calculate g20
num1=area_par0-bparu;
den1=area_par-bparu;

ratio1=num1/den1;

g20=ratio1;

%Calculate error in g20
n1=sqrt(area_par0+bpar_err^2);
d1=sqrt(area_par/npar+bpar_err^2);
r1=ratio1*sqrt((n1/num1)^2+(d1/den1)^2);
g20_err=r1;

end