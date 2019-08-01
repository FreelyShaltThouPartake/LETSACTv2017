function [AmsUtil, BmsUtil] = MSConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, MS, ConstrainedTech, option)
% Controls the growth and fall in sales of a given technology type

s = length(ConstrainedTech);
AmsUtil = sparse(s*analysis_horizon,N_Util*analysis_horizon*(Tcrit_Util+1));
BmsUtil = zeros(s*analysis_horizon,1);

switch option
    case 'NewSales'
        for k=1:analysis_horizon
            for i=1:s
                rowindex = i + (k-1)*s;
                AmsUtil(rowindex, squish(ConstrainedTech(i),k,1,N_Util,Tcrit_Util)) = ...
                    AmsUtil(rowindex, squish(ConstrainedTech(i),k,1,N_Util,Tcrit_Util)) - 1;
                for n=1:N_Util
                    AmsUtil(rowindex, squish(n,k,1,N_Util,Tcrit_Util)) = ...
                        AmsUtil(rowindex, squish(n,k,1,N_Util,Tcrit_Util)) + MS(ConstrainedTech(i),k);
                end
            end
        end
    case 'TotFleet'
        for k=1:analysis_horizon
            for i=1:s
                rowindex = i + (k-1)*s;
                AmsUtil(rowindex, squish(ConstrainedTech(i),k,1,N_Util,Tcrit_Util)) = ...
                    AmsUtil(rowindex, squish(ConstrainedTech(i),k,1,N_Util,Tcrit_Util)) - 1;
                for n=1:N_Util
                    AmsUtil(rowindex, squish(n,k,1,N_Util,Tcrit_Util)) = ...
                        AmsUtil(rowindex, squish(n,k,1,N_Util,Tcrit_Util)) + MS(ConstrainedTech(i),k);
                end
            end
        end
        for k=1:analysis_horizon
            for i=1:s
                rowindex = i + (k-1)*s;
                for j=1:Tcrit_Util
                    AmsUtil(rowindex,:) = AmsUtil(rowindex,:) - AstockUtil(squish(ConstrainedTech(i),j,k,N_Util,Tcrit_Util),:);
                    BmsUtil(rowindex) = BmsUtil(rowindex) - BstockUtil(squish(ConstrainedTech(i),j,k,N_Util,Tcrit_Util));
                end
                for n=1:N_Util
                    for j=1:Tcrit_Util
                        AmsUtil(rowindex,:) = AmsUtil(rowindex,:) + ...
                            MS(ConstrainedTech(i),k)*AstockUtil(squish(n,j,k,N_Util,Tcrit_Util),:);
                        BmsUtil(rowindex) = BmsUtil(rowindex) + ...
                            MS(ConstrainedTech(i),k)*BstockUtil(squish(n,j,k,N_Util,Tcrit_Util));
                    end
                end
            end
        end
end
end