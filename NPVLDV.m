function [fLDVnewyearly, fLDVoldyearly, fLDVretyearly, fLDVconstyearly, CF_LDV_total_new, CF_LDV_total_old]...
    = NPVLDV(N_LDV, Tcrit_LDV, analysis_horizon, FE_gas, FE_elec, FEdecline_rate, gasmiles, elecmiles, LDV_bodycosts, LDV_battcosts, LDV_maintcosts, LDV_retcosts, LDV_gascosts, LDV_eleccosts, AstockLDV, BstockLDV)
% Objective function that calculates net present value of all costs

%% Cost factor calculations
% Gas driving
[CF_LDV_new_gas, CF_LDV_old_gas, CF_LDV_ret] = createCFmatrix(N_LDV, Tcrit_LDV, analysis_horizon, FE_gas, FEdecline_rate, gasmiles, LDV_bodycosts, LDV_battcosts, LDV_maintcosts, LDV_retcosts, LDV_gascosts, 'gas');   % create 3-d matrix for cost factors; note that
                                                                     % costs for fuel use from electric driving are determined in the electric charging snippet  
                                                                     % All other LDV costs are accounted for in this instance of the createCFmatrix function
% Electric charging
[CF_LDV_new_elec, CF_LDV_old_elec, unused] = createCFmatrix(N_LDV, Tcrit_LDV, analysis_horizon, FE_elec, FEdecline_rate, elecmiles, 0, 0, 0, 0, LDV_eleccosts, 'elec');    % this invokation calculates fuel costs of electric driving only

CF_LDV_total_new = CF_LDV_new_gas + CF_LDV_new_elec;   % this cost factor for new LDVs includes deployment, maintenance, and fuel costs from gas+elec
CF_LDV_total_old = CF_LDV_old_gas + CF_LDV_old_elec;   % this cost factor for old LDVs includes deployment, maintenance, and fuel costs from gas+elec

%% Cost function setup
fLDVnewyearly = zeros(analysis_horizon,N_LDV*analysis_horizon*(Tcrit_LDV+1));
fLDVoldyearly = zeros(analysis_horizon,N_LDV*analysis_horizon*(Tcrit_LDV+1));
fLDVretyearly = zeros(analysis_horizon,N_LDV*analysis_horizon*(Tcrit_LDV+1));

fLDVconst = 0;
fLDVconstyearly = zeros(analysis_horizon,1);

% LDV emission constraints, market share constraints, vehicle #s and costs
newLDVend_index = N_LDV*analysis_horizon;
for k=1:analysis_horizon
    for i=1:N_LDV
        % Cost factors for new LDVs
        fLDVnewyearly(k,0 + squish(i,k,1,N_LDV,Tcrit_LDV)) = fLDVnewyearly(k,0 + squish(i,k,1,N_LDV,Tcrit_LDV)) + CF_LDV_total_new(i,k);
        for j=1:Tcrit_LDV         
            % Cost factors for old LDVs (note that (Ax - B)*CF will give
            % the total cost of old(i,j,k). However, optimal solution of
            % fLDV'x and fLDV'x + C are the same and thus C is ignored in the
            % formula for cost of old(i,j,k)
            fLDVoldyearly(k,:) = fLDVoldyearly(k,:) + AstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV),:)*CF_LDV_total_old(i,j,k);
            fLDVconst = fLDVconst + BstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV))*CF_LDV_total_old(i,j,k);
            fLDVconstyearly(k) = fLDVconstyearly(k) + BstockLDV(squish(i,j,k,N_LDV,Tcrit_LDV))*CF_LDV_total_old(i,j,k);    % store const cost values for yearly analysis later
            % Cost factors for retired LDVs
            fLDVretyearly(k,newLDVend_index + squish(i,j,k,N_LDV,Tcrit_LDV)) = fLDVretyearly(k,newLDVend_index + squish(i,j,k,N_LDV,Tcrit_LDV)) + CF_LDV_ret(i,j,k);
        end
    end
end

end

