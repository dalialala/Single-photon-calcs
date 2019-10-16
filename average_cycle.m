function [av_vec ,counter] = average_cycle(t, tn, t0, dt, pt, vec)
%This function finds the average pulse shape for a cycle of duration 
%pt with 2*tn cycles from a vector "vec"

%Input Variables
%   tn: pprox half number of cycles in time vector
%   t0: index of t=0
%   dt: timing resolution
%   pt: cycle duration in us
%   vec: vector with 2*tn+1 cycles

%Output variables
%   av_vec: average cycle
%   counter: number if times the cycle was averaged


%Number of bins in one cycle
t5u=round(pt/dt);

%Index where first cycle ends inmediately after t=0
tpt=t0+t5u-1;

%First cycle between t0 and tpt
av_vec=vec(t0:tpt);
%length(av_vec);

%Initialize counter
counter=1;

%Positive times
for i=1:tn-1
    av_vec=av_vec+vec(t0+i*t5u:tpt+i*t5u);
    t(t0+i*t5u);
    t(tpt+i*t5u);
    counter=counter+1;
end

%Negative times
for i=1:tn-1
    av_vec=av_vec+vec(t0-i*t5u+2:tpt-i*t5u+2);
    t(t0-i*t5u);
    t(tpt-i*t5u);
    counter=counter+1;
end
counter;
av_vec=av_vec/counter;
end

