function [Aclean, Bclean] = CleanTechConstraintSetup(N, T, Y, start_year, ref_index, dep_indices, ratios)
%Function to setup relationship between different forms of clean energy
%technologies capacity growth pegged to the growth of the reference technology; function
%setup as >= type

s = length(dep_indices);
Aclean = sparse(s*Y,N*Y*(T+1));
Bclean = zeros(s*Y,1);

for k=start_year:Y
    for i=1:s
        rowindex = i + (k-1)*s;
        Aclean(rowindex, squish(dep_indices(i),k,1,N,T)) = Aclean(rowindex, squish(dep_indices(i),k,1,N,T)) + 1;
        Aclean(rowindex, squish(ref_index,k,1,N,T)) = Aclean(rowindex, squish(ref_index,k,1,N,T)) - ratios(i);
    end
end
end