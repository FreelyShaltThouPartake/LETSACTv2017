function [AemLDV, BemLDV, EF_LDV_total_new, EF_LDV_total_old] = ...
    EmissionConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, gasmiles, elecmiles, FE_gas, FE_elec, FEdecline_rate, AstockLDV, BstockLDV, elecCO2content)
% Sets up emission constraint based on percentage reduction in emissions by
% 2050 relative to 1990 levels

%% Emission factor calculations
% LDV gas driving
gasCO2content = 8.92e-3;    % 8.92e-3 = tCO2/gal of gas
[EF_LDV_gas] = createEFmatrix(N_LDV, Tcrit_LDV, analysis_horizon, FE_gas, FEdecline_rate, gasmiles, gasCO2content,'gas');   % create 3-d matrix for emission factors
EF_LDV_gas_new = squeeze(EF_LDV_gas(:,1,:));    %EF for new vehicles separated from old vehicles
EF_LDV_gas_old = EF_LDV_gas(:,2:end,:);

% LDV electric charging
[EF_LDV_elec] = createEFmatrix(N_LDV, Tcrit_LDV, analysis_horizon, FE_elec, FEdecline_rate, elecmiles, elecCO2content,'elec');    % this invokation calculates emissions from electric driving only
EF_LDV_elec_new = squeeze(EF_LDV_elec(:,1,:));    %EF for new vehicles separated from old vehicles
EF_LDV_elec_old = EF_LDV_elec(:,2:end,:);
EF_LDV_total_old = EF_LDV_elec_old + EF_LDV_gas_old;  % combined EF from gas and electric driving
EF_LDV_total_new = EF_LDV_gas_new + EF_LDV_elec_new;

%% Constraint definition
AemLDV = sparse(analysis_horizon, N_LDV*analysis_horizon*(Tcrit_LDV+1));
% AemLDV = spalloc(Y, N_LDV*Y*(Tcrit_LDV+1),49280);
BemLDV = zeros(analysis_horizon,1);

for k=1:analysis_horizon
    for i=1:N_LDV
        % Emission factors for new LDVs
        AemLDV(k,0 + squish(i,k,1,N_LDV,Tcrit_LDV)) = AemLDV(k,0 + squish(i,k,1,N_LDV,Tcrit_LDV)) + EF_LDV_total_new(i,k); % squish function has third argument = 1 for new units
        for j=1:Tcrit_LDV
            % Emission factors for old LDVs
            AemLDV(k,:) = AemLDV(k,:) + AstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV),:)*EF_LDV_total_old(i,j,k);
            BemLDV(k) = BemLDV(k) + BstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV))*EF_LDV_total_old(i,j,k);
        end
    end
end
end

