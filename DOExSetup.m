function [LHC] = DOExSetup(param_num, sample_num)
% function sets up design of experiments using latin hypercubes
LHC = lhsdesign(sample_num,param_num)'; % transpose to get parameter values within each experiment in a given column
LHC(:,1) = 0.5;

[m,n] = size(LHC);

for j=1:n
    for i=1:m
        if LHC(i,j)<0.33
            LHC(i,j) = 0;
        elseif LHC(i,j)>0.66
            LHC(i,j) = 1;
        else
            LHC(i,j) = 0.5;
        end
    end
end
end

