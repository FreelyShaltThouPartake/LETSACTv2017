function [AdemLDV, BdemLDV] = DemandConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, AstockLDV, BstockLDV, option)
% Sets up demand constraint for vehicles miles traveled, expressed as
% number of vehicles since the miles/vehicle/year is assumed constant

AdemLDV = sparse(analysis_horizon, N_LDV*analysis_horizon*(Tcrit_LDV+1));
% AdemLDV = spalloc(Y, length(X0LDV),49280);
BdemLDV = zeros(analysis_horizon,1);

switch option
    case 'TotFleet'
        for k=1:analysis_horizon
            for i=1:N_LDV
                % Add new LDVs to total LDV count
                AdemLDV(k, squish(i,k,1,N_LDV,Tcrit_LDV)) = AdemLDV(k, squish(i,k,1,N_LDV,Tcrit_LDV)) + 1;
                for j=1:Tcrit_LDV
                    % Sum all old LDV numbers
                    AdemLDV(k,:) = AdemLDV(k,:) + AstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV),:);
                    BdemLDV(k) = BdemLDV(k) + BstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV));
                end
            end
        end
    case 'NewSales'
        for k=1:analysis_horizon
            for i=1:N_LDV
                % Add new LDVs to total LDV count
                AdemLDV(k, squish(i,k,1,N_LDV,Tcrit_LDV)) = AdemLDV(k, squish(i,k,1,N_LDV,Tcrit_LDV)) + 1;
            end
        end
end

