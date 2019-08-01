function [ARPS, BRPS] = RPSSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, RPS_value, RenewableTech, RPS_year, analysis_year0)
% Controls the total generation from renewable technologies as a fraction
% of total generation. Constraint set up as >=, use '-' sign in model to
% implement as standard <= type constraint

s = length(RenewableTech);
RPS_yearindex = RPS_year - analysis_year0;
years = analysis_horizon - RPS_yearindex + 1;
ARPS = sparse(years,N_Util*analysis_horizon*(Tcrit_Util+1));
BRPS = zeros(years,1);

for k=RPS_yearindex:analysis_horizon
    k
    for i=1:s
        ARPS(k - RPS_yearindex + 1, squish(RenewableTech(i),k,1,N_Util,Tcrit_Util)) = ...
            ARPS(k - RPS_yearindex + 1, squish(RenewableTech(i),k,1,N_Util,Tcrit_Util)) + 1;
        for j=1:Tcrit_Util
            ARPS(k - RPS_yearindex + 1,:) = ARPS(k - RPS_yearindex + 1,:) + AstockUtil(squish(RenewableTech(i),j,k,N_Util,Tcrit_Util),:);
            BRPS(k - RPS_yearindex + 1) = BRPS(k - RPS_yearindex + 1) + BstockUtil(squish(RenewableTech(i),j,k,N_Util,Tcrit_Util));
        end
    end
    for i=1:N_Util
        ARPS(k - RPS_yearindex + 1, squish(i,k,1,N_Util,Tcrit_Util)) = ARPS(k - RPS_yearindex + 1, squish(i,k,1,N_Util,Tcrit_Util)) - RPS_value;
        for j=1:Tcrit_Util
            ARPS(k - RPS_yearindex + 1,:) = ARPS(k - RPS_yearindex + 1,:) - RPS_value*AstockUtil(squish(i,j,k,N_Util,Tcrit_Util),:);
            BRPS(k - RPS_yearindex + 1) = BRPS(k - RPS_yearindex + 1) - RPS_value*BstockUtil(squish(i,j,k,N_Util,Tcrit_Util));
        end
    end
end
end