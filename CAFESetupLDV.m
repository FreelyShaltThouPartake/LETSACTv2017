function [Acafe, Bcafe] = CAFESetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, FE_combined, FE_cafe)
% Setup CAFE constraint

Acafe = sparse(analysis_horizon, N_LDV*analysis_horizon*(Tcrit_LDV+1));
% Acafe = spalloc(Y,length(X0LDV),49280);
Bcafe = zeros(analysis_horizon,1);

for k=1:analysis_horizon
    for i=1:N_LDV
        Acafe(k, 0 + squish(i,k,1,N_LDV,Tcrit_LDV)) = Acafe(k, 0 + squish(i,k,1,N_LDV,Tcrit_LDV)) + FE_combined(i,k) - FE_cafe(k);
        % in a given year k, sum(xi*FEi) >= (sum(xi))*CAFE
    end
    
end

end

