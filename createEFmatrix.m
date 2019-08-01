function [EFmatrix] = createEFmatrix(N, T, Y, FE, FEdecline_rate, miles, carb_intensity, opt)

EFmatrix = zeros(N,T+1,Y);    % Emission factor matrix (size in y dimension is T+1 because the 1st entry is for year = 0 i.e. new technologies that year)

if strcmp(opt,'util')==1
    HRnew = squeeze(FE(:,1,:)); % Heat rates of new plants by year of construction
    HRnew(1,:) = HRnew(1,:)*1.95;   % To account for coal CCS
    HRold = squeeze(FE(:,2:end,1)); % Heat rates of plants in 2014 initial fleet by age
end

for y=1:Y
    for tech=1:N
        if strcmp(opt,'gas')==1
            for age=1:T+1
                EFmatrix(tech,age,y) = FE(tech,26+y-age+1)*((1+FEdecline_rate)^age)*miles(tech,y)*carb_intensity;  % gpm and kWhpm increase by FEdecline_rate with age
            end
        elseif strcmp(opt,'elec')==1
            for age=1:T+1
                EFmatrix(tech,age,y) = FE(tech,26+y-age+1)*((1+FEdecline_rate)^age)*miles(tech,y)*carb_intensity(y);
            end
        elseif strcmp(opt,'util')==1
            if tech==1  % To account for coal CCS
                EFmatrix(tech,1,y) = HRnew(tech,y)*1.34*1e-3*carb_intensity(tech)*0.15*1e-3;    % To account for coal CCS
            else
                EFmatrix(tech,1,y) = HRnew(tech,y)*1e-3*carb_intensity(tech)*1e-3;    %First 1e-3 converts Btu/kWh into mmBtu/MWh, and the second converts kg CO2/mmBTu to tonnes/mmBtu
            end
            for age=2:T+1
                if age<=y
                    if tech==1
                        EFmatrix(tech,age,y) = HRnew(tech,y)*1.34*((1+FEdecline_rate)^age)*1e-3*carb_intensity(tech)*0.15*1e-3;   % heat rate increases by FEdecline_rate with age
                    else
                        EFmatrix(tech,age,y) = HRnew(tech,y)*((1+FEdecline_rate)^age)*1e-3*carb_intensity(tech)*1e-3;   % heat rate increases by FEdecline_rate with age
                    end
                else
                    EFmatrix(tech,age,y) = HRold(tech,age-y)*((1+FEdecline_rate)^age)*1e-3*carb_intensity(tech)*1e-3;
                end
            end
        end
    end
end
end