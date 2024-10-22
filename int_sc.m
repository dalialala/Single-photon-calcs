function [tvec, gate_corr, b_area, b_err] = int_sc(t, numer_g2, pt, tp, winc, file_list, HOM)

dt=(t(2)-t(1));
n=round(pt/dt);
t0=find(t>=0,1);
t=t(t0:t0+n-1);

igate=round(winc/dt);
%indext=76;
index=round((tp+winc)/dt)+1;
t(index);
indext=round((pt-t(index))/dt);
t(indext);

if winc ==0
    index=n;
    igate=n-1;
    indext=1;
end

tn=round(t(end)/pt);

singles3=double(numer_g2(2,:));
singles5=double(numer_g2(3,:));

[pulse3 ,counter] = average_cycle(t, tn, t0, dt, pt, singles3);
[pulse5 ,counter] = average_cycle(t, tn, t0, dt, pt, singles5);


nexp=length(file_list)*0.6/(pt*1e-6);

[~, ~, pba3, ~, pba5,  ~, pbce, pulse_re3, pulse_re5, ~, ~, ~, ~] = calc_backg(pulse3, pulse5, dt,  1, tp, winc, counter, nexp, HOM);
back=zeros(n,2*n);

% figure
% plot(fliplr(pulse_re5))
% 
% figure
% plot(pulse_re3)
for i=1:n
    for j=1:2*n
        if i<=index && i>=index-igate
            if j-i+1<=indext+igate && j-i+1>=indext
                back(i,j)=pba5;
            end
        else
            pulse_re3(i)=0;
            pulse_re5(i)=0;
        end
    end
end

b3b5=back*pba3;

if winc ==0
    back=pba5*ones(n,2*n);
    b3b5=back*pba3;
    index=2*n;
    igate=2*n-1;
end


p3b5=pulse_re3'.*back;

pulse5m=zeros(n,2*n);
vec0=zeros(1,n);
pulse5m(1,:)=[fliplr(pulse_re5), vec0];

for i=2:n
    if i<=index && i>=index-igate
        pulse5m(i, :)=[vec0(1:i-1), fliplr(pulse_re5), vec0(i:end)];
    end
end

p5b3=pulse5m*pba3;

p3p5=((pulse_re3))'.*pulse5m;

tot=(p3b5+p5b3+b3b5)*nexp;

totp=(p3b5+p5b3+b3b5+p3p5)*nexp;

% [T1, Tau]=meshgrid([t, t+t(end)], t);
% figure
% surf(T1, Tau,totp), colorbar, view(2);
% caxis([0 5])
% shading interp


tvec=-n*dt:dt:(n-1)*dt;
[T1, Tau]=meshgrid(tvec, t);
figure
surf(T1, Tau,tot),view(2);
shading interp
% xlim([-3, 3])
% ylim([2.5, 5.0])
caxis([0 5])
pbaspect([1 1 1])

gate_corr=zeros(1, 2*n);
dummy=sum(tot,1);

if winc ==0
    gate_corr=dummy;
else
    gate_corr(n-igate:n+igate)=dummy(n-igate:n+igate);
end

b_area=sum(gate_corr);
b_err=pbce*(2*igate+1);

% figure
% plot(tvec, gate_corr)
end