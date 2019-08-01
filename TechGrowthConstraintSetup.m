function [AnewLDV, BnewLDV] = TechGrowthConstraintSetup(N_LDV, Tcrit_LDV, analysis_horizon, MS, initcap)
% Controls the growth and fall in sales of a given technology type

pnew_inc = MS(1:N_LDV,:);
pnew_dec = MS(N_LDV+1:end,:);

apply_inc_const_to_tech = [2,3,4];
apply_dec_const_to_tech = 1;

AnewLDV = sparse(2*N_LDV*analysis_horizon,N_LDV*analysis_horizon*(Tcrit_LDV+1));
% AnewtotbytechLDV = spalloc(2*N_LDV*Y,length(X0LDV),961);
BnewLDV = zeros(2*N_LDV*analysis_horizon,1);

for k=1:analysis_horizon
    for i=1:length(apply_inc_const_to_tech)
        tech = apply_inc_const_to_tech(i);
        rowindex = tech + (k-1)*N_LDV;
        if k==1
            AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) = AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) + 1;
            BnewLDV(rowindex) = BnewLDV(rowindex) + (1 + pnew_inc(tech,k))*initcap(tech);
        else
            AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) = AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) + 1;
            AnewLDV(rowindex, squish(tech,k-1,1,N_LDV,Tcrit_LDV)) = AnewLDV(rowindex, squish(tech,k-1,1,N_LDV,Tcrit_LDV)) - (1 + pnew_inc(tech,k));
        end
    end
end

for k=1:analysis_horizon
    for i=1:length(apply_dec_const_to_tech)
        tech = apply_dec_const_to_tech(i);
        rowindex = N_LDV*analysis_horizon + tech + (k-1)*N_LDV;
        if k==1
            AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) = AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) - 1;
            BnewLDV(rowindex) = BnewLDV(rowindex) + (1 - pnew_dec(tech,k))*initcap(tech);
        else
            AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) = AnewLDV(rowindex, squish(tech,k,1,N_LDV,Tcrit_LDV)) - 1;
            AnewLDV(rowindex, squish(tech,k-1,1,N_LDV,Tcrit_LDV)) = AnewLDV(rowindex, squish(tech,k-1,1,N_LDV,Tcrit_LDV)) + (1 - pnew_dec(tech,k));
        end
    end
end

end

