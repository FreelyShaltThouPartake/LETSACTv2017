function [fUtilnewyearly, fUtiloldyearly, fUtilretyearly, fUtilconstyearly, CFnew_Util, CFold_Util, CFret_Util] = NPVUtil(N_Util, Tcrit_Util, analysis_horizon, Util_heatrates, FEdecline_rate,...
    Util_techcosts, Util_maintcosts, Util_retcosts, Util_fuelprices, AstockUtil, BstockUtil)

%% Cost factor calculations

[CFnew_Util, CFold_Util, CFret_Util] = createCFmatrix(N_Util, Tcrit_Util, analysis_horizon, Util_heatrates, FEdecline_rate, 0, Util_techcosts, 0, Util_maintcosts, Util_retcosts, Util_fuelprices, 'util');

fUtilconst = 0;
fUtilconstyearly = zeros(analysis_horizon,1);
fUtilnewyearly = zeros(analysis_horizon,N_Util*analysis_horizon*(1+Tcrit_Util));
fUtiloldyearly = zeros(analysis_horizon,N_Util*analysis_horizon*(1+Tcrit_Util));
fUtilretyearly = zeros(analysis_horizon,N_Util*analysis_horizon*(1+Tcrit_Util));

for k=1:analysis_horizon
    for i=1:N_Util
        % Cost factors for new Utils
        fUtilnewyearly(k,0 + squish(i,k,1,N_Util,Tcrit_Util)) = fUtilnewyearly(k,0 + squish(i,k,1,N_Util,Tcrit_Util)) + CFnew_Util(i,k);
        for j=1:Tcrit_Util
            % Cost factors for old Utils (note that (Ax - C)*CFUtil will give
            % the total cost of old(i,j,k). However, optimal solution of
            % fUtil'x and fUtil'x + C are the same and thus C is ignored in the
            % formula for cost ofUtil old(i,j,k)
            fUtiloldyearly(k,:) = fUtiloldyearly(k,:) + AstockUtil(squish(i,j,k,N_Util,Tcrit_Util),:)*CFold_Util(i,j,k);
            fUtilconst = fUtilconst + BstockUtil(squish(i,j,k,N_Util,Tcrit_Util))*CFold_Util(i,j,k);
            fUtilconstyearly(k) = fUtilconstyearly(k) + BstockUtil(squish(i,j,k,N_Util,Tcrit_Util))*CFold_Util(i,j,k);    % store const cost values for yearly analysis later

            % Cost factors for retired Utils
            fUtilretyearly(k,N_Util*analysis_horizon + squish(i,j,k,N_Util,Tcrit_Util)) = fUtilretyearly(k,N_Util*analysis_horizon + squish(i,j,k,N_Util,Tcrit_Util)) + CFret_Util(i,j,k);
        end
    end
end
end