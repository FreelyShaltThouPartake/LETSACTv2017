function [CnewMatrix, ColdMatrix, CretMatrix] = createCFmatrix(N, T, Y, FE, FEdecline_rate, miles, C_body, C_batt, C_maint, C_ret, FuelPrice, opt)

CnewMatrix = zeros(N,Y);
ColdMatrix = zeros(N,T,Y);
CretMatrix = zeros(N,T,Y);

if strcmp(opt,'util')==1
    % FE here is plant heat rate in Btu/kWh, and increase in heat rate
    % is built into the spreadsheet from which data is read
    % Fuel price is in $/mmBtu
    HRnew = squeeze(FE(:,1,:)); % Heat rates of new plants by year of construction
    HRold = squeeze(FE(:,2:end,1)); % Heat rates of plants in 2014 initial fleet by age
%     fin_rate = 0.2;   % Power plant financing rate over 20 years including ROI (12%), risk (3%), and tax (5%)
%     fin_period = 20;
%     annuity_factor = ((1+fin_rate)^fin_period - 1)/(fin_rate*(1 + fin_rate)^fin_period);
end

for y=1:Y
    if strcmp(opt,'gas')==1
        % FE here is fuel economy in mpg
        % Fuel price is in $/gal
        for tech=1:N
            CnewMatrix(tech,y) = C_body(tech,y) + C_batt(tech,y) + miles(tech,y)*FuelPrice(tech,y)*FE(tech,25+y+1) + C_maint(tech,y);
            for age=1:T
                ColdMatrix(tech,age,y) = miles(tech,y)*FuelPrice(tech,y)*FE(tech,25+y-age+1)*((1+FEdecline_rate)^age) + C_maint(tech,y);
                % first term in Cold Matrix is fuel cost and second is maintenance cost
                % FE for LDVs is assumed to reduce by 0.5% with age in any given year
                
                CretMatrix(tech,age,y) = C_ret(tech,age);
            end
        end
        
    elseif strcmp(opt,'elec')==1
        for tech=1:N
            CnewMatrix(tech,y) = FE(tech,24+y+1)*miles(tech,y)*FuelPrice(tech,y);
            for age=1:T
                ColdMatrix(tech,age,y) = FE(tech,2014+abs(y-age)-1989+1)*((1+FEdecline_rate)^age)*miles(tech,y)*FuelPrice(tech,y);
                % only fuel cost is needed for electric driving since other
                % cost factors are obtained when this function is run for
                % gasoline
                
                CretMatrix = 0;
            end
        end
        
    elseif strcmp(opt,'util')==1
        for tech=1:N
            if tech==1  % To account for CCS on new PC plants
                CnewMatrix(tech,y) = C_body(tech,y)*1.95 + FuelPrice(tech,y)*HRnew(tech,y)*1e-3*1.34 + C_maint(tech,y)*1.63;
            else
                CnewMatrix(tech,y) = C_body(tech,y) + FuelPrice(tech,y)*HRnew(tech,y)*1e-3 + C_maint(tech,y);
            end
            % Capital costs assumed to be overnight, and implemented as
            % the upfront total payment paid over the first year
            % Conversion factors in the middle term are for getting $/MWh
            % from $/mmBtu and Btu/kWh
            for age=1:T
                if age<=y
                    if tech==1
                        ColdMatrix(tech,age,y) = FuelPrice(tech,y)*HRnew(tech,y)*1.34*((1+FEdecline_rate)^age)*1e-3 + C_maint(tech,y)*1.63;
                        % first term in Cold Matrix is fuel cost and second is maintenance cost
                    else
                        ColdMatrix(tech,age,y) = FuelPrice(tech,y)*HRnew(tech,y)*((1+FEdecline_rate)^age)*1e-3 + C_maint(tech,y);
                        % first term in Cold Matrix is fuel cost and second is maintenance cost
                    end
                    CretMatrix(tech,age,y) = C_ret(tech,age,y);
                else
                    ColdMatrix(tech,age,y) = (FuelPrice(tech,y)*HRold(tech,age-y)*((1+FEdecline_rate)^age)*1e-3 + C_maint(tech,y));
                    % first term in Cold Matrix is fuel cost and second is maintenance cost
                    CretMatrix(tech,age,y) = C_ret(tech,age,y);
                end
            end
        end
    end
end
