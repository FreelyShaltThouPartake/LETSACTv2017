function [AemUtil, BemUtil, EFnew_Util, EFold_Util] = EmissionConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon,...
    HeatRates, FEdecline_rate, carb_intensity, AstockUtil, BstockUtil)
    
[EFUtil] = createEFmatrix(N_Util, Tcrit_Util, analysis_horizon, HeatRates, FEdecline_rate, 0, carb_intensity,'util');
EFnew_Util = squeeze(EFUtil(:,1,:));  % squeeze out emission factors for newly deployed units in all years
EFold_Util = EFUtil(:,2:end,:);   % emission factors for old units

AemUtil = sparse(analysis_horizon, N_Util*analysis_horizon*(1+Tcrit_Util));
BemUtil = zeros(analysis_horizon,1);
for k=1:analysis_horizon
    for i=1:N_Util
        % Emission factors for new Utils
        AemUtil(k,0 + squish(i,k,1,N_Util,Tcrit_Util)) = AemUtil(k,0 + squish(i,k,1,N_Util,Tcrit_Util)) + EFnew_Util(i,k);
        for j=1:Tcrit_Util         
            % Emission factors for old Utils
            AemUtil(k,:) = AemUtil(k,:) + AstockUtil(squish(i,j,k,N_Util,Tcrit_Util),:)*EFold_Util(i,j,k);
            BemUtil(k) = BemUtil(k) + BstockUtil(squish(i,j,k,N_Util,Tcrit_Util))*EFold_Util(i,j,k);
        end
    end
end

end