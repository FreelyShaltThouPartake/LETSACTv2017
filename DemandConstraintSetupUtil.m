function [AdemUtil, BdemUtil] = DemandConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, option)

AdemUtil = sparse(analysis_horizon, N_Util*analysis_horizon*(1+Tcrit_Util));
BdemUtil = zeros(analysis_horizon,1);

switch option
    case 'TotFleet'
        for k=1:analysis_horizon
            k
            for i=1:N_Util
                % Add new Utils to total Util count
                AdemUtil(k, squish(i,k,1,N_Util,Tcrit_Util)) = AdemUtil(k, squish(i,k,1,N_Util,Tcrit_Util)) + 1;
                for j=1:Tcrit_Util
                    % Sum all old Util numbers
                    AdemUtil(k,:) = AdemUtil(k,:) + AstockUtil(squish(i,j,k,N_Util,Tcrit_Util),:);
                    BdemUtil(k) = BdemUtil(k) + BstockUtil(squish(i,j,k,N_Util,Tcrit_Util));
                end
            end
        end
    case 'New'
        for k=1:analysis_horizon
            k
            for i=1:N_Util
                % Add new generation to total total generation
                AdemUtil(k, squish(i,k,1,N_Util,Tcrit_Util)) = AdemUtil(k, squish(i,k,1,N_Util,Tcrit_Util)) + 1;
            end
        end     
end

end