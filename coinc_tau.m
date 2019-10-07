function [mat_out] = coinc_tau(t1, matrix)
%This function shifts the coincidence space t1-t2 to t1-tau
%where tau=t1-t2

%Input Variables
%   matrix: P(t1,t2) coincidence matrix

%Output Variables
%   matrix: P(t1,tau) coincidence matrix

[n,m]=size(matrix);

mat_out=zeros(n, m);

tm=max(t1);

for i=1:n
    ti=t1(i);
    for j=1:m
        tj=t1(j);
        tau=(ti-tj);
        if tau<0
            tau=tm+tau;
        end
        indexj=find(t1 >=tau, 1);
        mat_out(i,j)=matrix(i, indexj);
    end
end

% for i=1:n
%     for j=1:n-1
%         mat_out(i,j)=matrix(i, j+1);
%     end
% end

end

