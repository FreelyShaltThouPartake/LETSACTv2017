function [AprodcapLDV, BprodcapLDV] = ProdCapConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon_LP, analysis_year0, climate_action_year, initprod)
% Restricts growth of total production capacity of vehicles per year

AprodcapLDV = sparse(analysis_horizon_LP, N_LDV*analysis_horizon_LP*(Tcrit_LDV+1));
BprodcapLDV = zeros(analysis_horizon_LP,1);

for k=1:analysis_horizon_LP
    for i=1:N_LDV
        AprodcapLDV(k, squish(i,k,1,N_LDV,Tcrit_LDV)) = AprodcapLDV(k, squish(i,k,1,N_LDV,Tcrit_LDV)) + 1;
    end
end

for k=1:analysis_horizon_LP
    BprodcapLDV(k) =  initprod*(1.005)^k;
end

for k=(climate_action_year - analysis_year0):analysis_horizon_LP
    if k<=(climate_action_year - analysis_year0) + 10
        BprodcapLDV(k) =  BprodcapLDV(k) + 1e6*k;
    else
        BprodcapLDV(k) =  BprodcapLDV(k) + 1e6*((climate_action_year - analysis_year0) + 10);
    end
end
end