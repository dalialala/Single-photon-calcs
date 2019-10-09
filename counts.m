function [corr, single, t1gate_vec, t1tgate_vec] = counts(tau_2, tau_3, numer_g3, tp, t1gate, pt)
%Convert g3 data to double
g3=squeeze(double(numer_g3));

%Size of g3
[m,n]=size(g3);


%Time vectors in us
t1=tau_3*1e6;
t=tau_2*1e6;

%Bin size
dt=t(2)-t(1);

%Single counts vector
single=sum(g3(:,1:end), 2);

%Spad Correlation vector
corr=sum(g3(1:end,:), 1);

%Time pulse starts in clock time
tlp = find(t1>=tp,1);
t1(tlp);

%Maximum pulse in clock time
tp_max = find(t1>=3.37,1);

%Time pulse finishes in clock time
%thp = find(t1>=3.47,1);
thp =find(t1>=tp+t1gate,1);
t1(thp);
thp-tlp+1;

%g3c=g3;
g3c=zeros(thp-tlp+1,n);
t1i=0;
t1f=0;

tmin = t(1);
tmax = t(end);
t0=find(t>=0,1);

tn=round(tmax/pt);

%Time that gate starts for each t1
tb=-3.4941+dt;
%Time that gate finishes for each t1 original 4.14
te=tb+t1gate;

(te)-(tb)

%For loop runs first on time t1
for i=tlp:thp
    %Second for loop runs for gating in tau
    for j=1:tn
        t1i=find(t-((j-1)*(pt)+t1(i)+tb)<=dt/2, 1, 'last');
        ti=t(t1i);
        t1f=find(t-((j-1)*(pt)+t1(i)+te)<=dt/2, 1, 'last')+1;
        t2=t(t1f);
        t1f-t1i+1;
        g3c(i-tlp+1, t1i:t1f)=g3(i,t1i:t1f);
    end
end


%Negative times
%for i=167:167
%    for j=1:2
for i=tlp:thp
    for j=1:tn-1
        t1i=find(t-(-j*pt+t1(i)+tb)<=dt/2, 1, 'last');
        ti=t(t1i);
        t1f=find(t-(-j*pt+t1(i)+te)<=dt/2, 1, 'last');
        t2=t(t1f);
        dtau=t2-ti;
        g3c(i-tlp+1, t1i:t1f)=g3(i,t1i:t1f);
    end
end

t1gate_vec=sum(g3(tlp:thp,:),1);
t1tgate_vec=sum(g3c(1:end,:),1);

%Plot
% figure
% semilogy(t, corr)
% hold
% semilogy(t, t1gate_vec)
% semilogy(t, t1tgate_vec)
% hold
% xlim([-10 10])
end
