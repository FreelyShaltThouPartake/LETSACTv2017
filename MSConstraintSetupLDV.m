function [AmsLDV, BmsLDV] = MSConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, AstockLDV, BstockLDV, MS, ConstrainedTech, option)
% Controls the growth and fall in sales of a given technology type

s = length(ConstrainedTech);
AmsLDV = sparse(s*analysis_horizon,N_LDV*analysis_horizon*(Tcrit_LDV+1));
BmsLDV = zeros(s*analysis_horizon,1);

switch option
    case 'NewSales'
        for k=1:analysis_horizon
            for i=1:s
                rowindex = i + (k-1)*s;
                AmsLDV(rowindex, squish(ConstrainedTech(i),k,1,N_LDV,Tcrit_LDV)) = ...
                    AmsLDV(rowindex, squish(ConstrainedTech(i),k,1,N_LDV,Tcrit_LDV)) - 1;
                for n=1:N_LDV
                    AmsLDV(rowindex, squish(n,k,1,N_LDV,Tcrit_LDV)) = ...
                        AmsLDV(rowindex, squish(n,k,1,N_LDV,Tcrit_LDV)) + MS(ConstrainedTech(i),k);
                end
            end
        end
    case 'TotFleet'
        for k=1:analysis_horizon
            for i=1:s
                rowindex = i + (k-1)*s;
                AmsLDV(rowindex, squish(ConstrainedTech(i),k,1,N_LDV,Tcrit_LDV)) = ...
                    AmsLDV(rowindex, squish(ConstrainedTech(i),k,1,N_LDV,Tcrit_LDV)) - 1;
                for n=1:N_LDV
                    AmsLDV(rowindex, squish(n,k,1,N_LDV,Tcrit_LDV)) = ...
                        AmsLDV(rowindex, squish(n,k,1,N_LDV,Tcrit_LDV)) + MS(ConstrainedTech(i),k);
                end
            end
        end
        for k=1:analysis_horizon
            for i=1:s
                rowindex = i + (k-1)*s;
                for j=1:Tcrit_LDV
                    AmsLDV(rowindex,:) = AmsLDV(rowindex,:) - AstockLDV(squish(ConstrainedTech(i),j,k,N_LDV,Tcrit_LDV),:);
                    BmsLDV(rowindex) = BmsLDV(rowindex) - BstockLDV(squish(ConstrainedTech(i),j,k,N_LDV,Tcrit_LDV));
                end
                for n=1:N_LDV
                    for j=1:Tcrit_LDV
                        AmsLDV(rowindex,:) = AmsLDV(rowindex,:) + ...
                            MS(ConstrainedTech(i),k)*AstockLDV(squish(n,j,k,N_LDV,Tcrit_LDV),:);
                        BmsLDV(rowindex) = BmsLDV(rowindex) + ...
                            MS(ConstrainedTech(i),k)*BstockLDV(squish(n,j,k,N_LDV,Tcrit_LDV));
                    end
                end
            end
        end
end
end