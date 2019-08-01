% (c) Sarang Supekar
% University of Michigan, Ann Arbor, MI, USA

tic
clc
clear all
close all
% warning('off','all');

% % initialize Java for xlwrite function (substitute for xlswrite, which does
% % not work on Mac)
% javaaddpath('poi_library/poi-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
% javaaddpath('poi_library/xmlbeans-2.3.0.jar');
% javaaddpath('poi_library/dom4j-1.6.1.jar');
% javaaddpath('poi_library/stax-api-1.0.1.jar');

%% Define basic model variables and parameters
disp('Defining basic model variables and parameters.');

N_LDV = 4;  % no. of vehicle technologies considered
Tcrit_LDV = 25; % Critical age after which vehicle is assumed to be scrapped
discountrate = 0.07;

reference_year = 2010;
target_year = 2050;
reduction_target = 0.71;    % %reduction in CO2 by target year
analysis_year0 = 2014;
analysis_yearend = target_year + Tcrit_LDV;
analysis_horizon = analysis_yearend - analysis_year0;  % no. of years considered in LP analysis
planning_horizon = target_year - (analysis_year0 + 1) + 1;
caylength = 4;

include_SCC = 0;   % binary variable for whether or not to include social cost of carbon in the objective function; 1 = include, 0 = exclude

%% Create design of experiments
disp('Creating design of experiments.');

param_num = 3;
levels = 3;
exp_num = param_num^levels;
LDVExps = fullfact(3*ones(1,param_num))' - 1;
% save('LDV_DOEx.mat','LDVExps');
CaseNames = {'Low','Nom','High'};
CaseNumbers = [0, 1, 2];
input_filename = 'LDV_Master_Input_File.xlsx';

X = zeros(N_LDV*analysis_horizon*(1 + Tcrit_LDV),caylength,exp_num);
Xbau = zeros(N_LDV*analysis_horizon*(1 + Tcrit_LDV),exp_num);
R = zeros(caylength,exp_num);
S = zeros(caylength,exp_num);
CACTION = zeros(caylength,exp_num);
CBAU = zeros(1,exp_num);
CINCR = zeros(caylength,exp_num);
CEFFECT = zeros(caylength,exp_num);
NEWTOT = zeros(N_LDV,caylength,exp_num);
OLDTOT = zeros(N_LDV,caylength,exp_num);
RETTOT = zeros(N_LDV,caylength,exp_num);
DEADTOT = zeros(N_LDV,caylength,exp_num);
NEWTOTBAU = zeros(N_LDV,exp_num);
OLDTOTBAU = zeros(N_LDV,exp_num);
RETTOTBAU = zeros(N_LDV,exp_num);
DEADTOTBAU = zeros(N_LDV,exp_num);
ECA = zeros(target_year - analysis_year0,caylength,exp_num);
EBAU = zeros(target_year - analysis_year0,exp_num);
EVCHARGEDEM = zeros(analysis_horizon,caylength,exp_num);
EVCHARGEDEMBAU = zeros(analysis_horizon,exp_num);

%% Setup vehicle stocks
disp('Setting up vehicle stocks.');
% Initial (2014) fleet composition
initfleetLDV = xlsread(input_filename,'Nom','B3:Z6'); % intial (2014) fleet of LDVs
inittotalbytypeLDV = sum(initfleetLDV,2);
inittotalLDV = sum(sum(initfleetLDV));
initnewLDV = initfleetLDV(:,1); % first row stores values for new vehicles in 2014
initoldLDV = initfleetLDV(:,2:end);   % old vehicles in 2010

% Vehicle survival probabilities
disc_prob_LDV = xlsread(input_filename,'Nom','B10:AA13'); % discard probabilities for LDVs of different ages
P_LDV = 1.-disc_prob_LDV;   % survival probabilities for LDVs of different ages

[AstockLDV, BstockLDV] = StockSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, initoldLDV, initnewLDV, P_LDV);
save('.\MATS\LDV\LDVMatsCommon.mat','AstockLDV','BstockLDV');

%% Annual demand constraint on vehicles sales or total vehicle stock (actual demand values need to be added to RHS based on scenario)
disp('Setting up demand constraint.');
[AdemLDV, BdemLDV] = DemandConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, AstockLDV, BstockLDV, 'TotFleet');
save('.\MATS\LDV\LDVMatsCommon.mat','AdemLDV','BdemLDV','-append');

%% Market share constraint
disp('Setting up market shares constraint.');
MS = xlsread(input_filename,'Nom','B96:BJ99');  % Projected market shares of new vehicle sales for BAU scenario
ConstrainedTech = [2,4];  % Index of technologies to apply >= style market share constraints
[AmsLDV, BmsLDV] = MSConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, AstockLDV, BstockLDV, MS, ConstrainedTech, 'NewSales'); % Do not use "-" sign before this >= type contraint in the optim function call since it is formulated as <=0 type within this setup function
save('.\MATS\LDV\LDVMatsCommon.mat','AmsLDV','BmsLDV','-append');

%% Emission constraint setup data
disp('Setting up data for emission constraint.');
reduction_percent_2050 = reduction_target;
emission_history = xlsread('EmissionConstraints.xls','EmissionHistoryLDV');
cuml_emission_lockin = 0;   % this variable is to be used to account for cumulative emission lock-ins as estimated by the SD model
em_2010 = emission_history(1); % CO2 emissions (tCO2) in 2010

% Ideal linear emission reduction curve equation
target_year_emissions = (1 - reduction_percent_2050)*em_2010;
ideal_em_curve_slope = (target_year_emissions - emission_history(1))/(target_year - reference_year);

linear_emcon = zeros(analysis_horizon + analysis_year0 - reference_year + 1,1);
for year=1:length(linear_emcon)
    if year<=target_year - reference_year + 1
        linear_emcon(year) = ideal_em_curve_slope*(year) + emission_history(1);
    else
        linear_emcon(year) = ideal_em_curve_slope*(target_year - reference_year) + emission_history(1);
    end
end

% Calculate emission budget used up from reference year to year 0
history_years = reference_year:analysis_year0;  % time array from 2010 to 2014
emission_debt = sum(emission_history - (ideal_em_curve_slope.*(history_years - reference_year) + emission_history(1))'); % Sum of emission deficits from reference year to initial year assuming reference case as climate action in reference year followed by a linear emission reduction till 2050
history_years = analysis_year0+1:target_year;
gross_emission_budget = sum(ideal_em_curve_slope.*(history_years - reference_year) + emission_history(1)); % sum allowable emission levels from analysis year 1 to target year
net_emission_budget = gross_emission_budget - emission_debt;

for run=14:14
    fprintf('\n\nRunning experiment # %d...\n\n',run)
    %% create output and mat filenames and write input data parameter values for experiment # run
    excel_filename = strcat('\Results\LDV\CAD_LDV_',num2str(run),'.xlsx');
    mat_filename = strcat('.\MATS\LDV\CAD_LDV_Mats_',num2str(run),'.mat');
    
    %% List of parameters for uncertainty analysis
    % 
    % 1	Energy Prices, vehicle cap and ret costs, PHEV:BEV market share
    % 2	Fuel economy, decrease in fuel economy with age, electric charging
    % emissions
    % 3	Projected vehicle population, miles/year, ratio of electric to
    % non-electric driving
    
    %% Load essential data from Excel files
    disp('Loading essential data from Excel files.');

    % Vehicle costs and fuel economy
    FE_gas = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(2,run))),'B50:CJ53'); % Fuel economy of gasoline driving in gal/mile as a function of time
    FE_elec = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(2,run))),'B57:CJ60');    % Fuel economy of electric charging from grid in MWh/mile as a function of time
    FE_combined = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(2,run))),'AB64:CJ67');    % Combined fuel economy in MPGe as a function of time
    FE_cafe = xlsread(input_filename,'Nom','AB71:CJ71');    % CAFE standards in mpg
    
    LDV_bodycosts = xlsread(input_filename,'Nom','B22:BJ25');   % N_LDV x Y matrix containing vehicle body costs without battery ($/vehicle)
    LDV_battcosts = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(1,run))),'B29:BJ32');   % N_LDV x Y matrix containing vehicle battery costs only ($/vehicle)
    LDV_maintcosts = xlsread(input_filename,'Nom','B36:BJ39');   % N_LDV x Y vector containing maintenance costs ($/vehicle/year)
    LDV_retcosts = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(1,run))),'B43:AA46'); % N_LDV x Tcrit_LDV vector containing retirement costs (or vehicle value) as a function of age ($/vehicle/year)
    LDV_gascoststemp = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(1,run))),'B17:CS17');  
    LDV_gascosts = [LDV_gascoststemp;LDV_gascoststemp;LDV_gascoststemp;LDV_gascoststemp];   % N_LDV x Y matrix containing fossil fuel price trends for all techs
    LDV_eleccoststemp = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(1,run))),'B18:CS18')*10; % N_LDV x Y matrix containing electric charging price trends for all techs
    LDV_eleccosts = [LDV_eleccoststemp;LDV_eleccoststemp;LDV_eleccoststemp;LDV_eleccoststemp];
    FEdecline_rate = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(2,run))),'B117');  % Annual rate at which gpm and MWhpm increase with every year of vehicle ageing
%     SCC = xlsread(input_filename,'SCC');   %Social cost of carbon in $/ton of CO2
    SCC = 100;
    
    totfleetdemand = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(3,run))),'B103:BJ103')'; % Growth in vehicle population over analysis time horizon    
    miles = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(3,run))),'B75:BJ78');
    elec_driving_vector = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(3,run))),'B82:B85');   %   divide miles into electric driving and gas driving miles
    for k=1:analysis_horizon
        for i=1:N_LDV
            elecmiles(i,k) = elec_driving_vector(i)*miles(i,k);
        end
    end
    gasmiles = miles - elecmiles;
    
    % N_LDV x Y matrix containing daily EV charging load in kWh
    Daily_Charging_Load = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(3,run))),'B111:BJ114');
    save(mat_filename,'gasmiles','elecmiles');
   
    %% Retirement constraint
    % disp('Setting up retirement constraint.');
    % min_retirement_age = 11;
    %
    % if min_retirement_age>1
    %     for k=1:analysis_horizon
    %         for i=1:N_LDV
    %             for j=1:min_retirement_age-1
    %                 ubLDV(N_LDV*analysis_horizon + squish(i,j,k,N_LDV,Tcrit_LDV)) = 0;
    %             end
    %         end
    %     end
    % end
    
    %% NPV objective function
    disp('Setting up objective function.');
    
    [fLDVnewyearly, fLDVoldyearly, fLDVretyearly, fLDVconstyearly, CF_LDV_total_new, CF_LDV_total_old] = NPVLDV(N_LDV, Tcrit_LDV, analysis_horizon, FE_gas, FE_elec, FEdecline_rate, gasmiles, elecmiles,...
        LDV_bodycosts, LDV_battcosts, LDV_maintcosts, LDV_retcosts, LDV_gascosts, LDV_eleccosts, AstockLDV, BstockLDV);
    
    fLDV = zeros(1,N_LDV*analysis_horizon*(Tcrit_LDV+1));
    for k=1:analysis_horizon
        fLDVnewyearly(k,:) = fLDVnewyearly(k,:)/(1 + discountrate)^(k-1);
        fLDVoldyearly(k,:) = fLDVoldyearly(k,:)/(1 + discountrate)^(k-1);
        fLDVretyearly(k,:) = fLDVretyearly(k,:)/(1 + discountrate)^(k-1);
        fLDVconstyearly(k,:) = fLDVconstyearly(k,:)/(1 + discountrate)^(k-1);
        fLDV(1,:) = fLDV(1,:) + fLDVnewyearly(k,:) + fLDVoldyearly(k,:) + fLDVretyearly(k,:);
    end
    save(mat_filename, 'fLDV','fLDVnewyearly', 'fLDVoldyearly', 'fLDVretyearly', 'fLDVconstyearly', 'CF_LDV_total_new', 'CF_LDV_total_old','-append');
    
    %% CAFE constraint
    disp('Setting up CAFE constraint.');
    [Acafe, Bcafe] = CAFESetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, FE_combined, FE_cafe);
    save(mat_filename,'Acafe','Bcafe','-append');
    
    %% Create constraint for proportionate BEV deployment
    disp('Creating proportional EV deployment constraint');
    EVratio = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(3,run))),'B89:B92');
    ref_index = 4;  % PHEV is used as reference
    dep_indices = [3];    % BEV used as dependent
    ratios = EVratio(dep_indices);    %
    [AcleanLDV, BcleanLDV] = CleanTechConstraintSetup(N_LDV, Tcrit_LDV, analysis_horizon, 2015 - analysis_year0, ref_index, dep_indices, ratios);
    save(mat_filename,'AcleanLDV','BcleanLDV','-append');
    
    %% BAU emission factors setup
    disp('Setting up BAU emission factors.');
    elecCO2contentbau = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(2,run))),'B107:BJ107'); % Capacity weighted average emission intensity tCO2/MWh as a function of time and CAY using an iterative process with Util model
    elecCO2contentmat = xlsread(input_filename,char(CaseNames(CaseNumbers==LDVExps(2,run))),'B121:BJ135'); % Capacity weighted average emission intensity tCO2/MWh as a function of time and CAY using an iterative process with Util model
    
    % Compute A matrix and B vector for emission constraint
    [AemLDVbau, BemLDVbau, EF_LDV_total_new_bau, EF_LDV_total_old_bau] =...
        EmissionConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, gasmiles, elecmiles, FE_gas, FE_elec, FEdecline_rate, AstockLDV, BstockLDV, elecCO2contentbau);
    save(mat_filename,'AemLDVbau','BemLDVbau','EF_LDV_total_new_bau','EF_LDV_total_old_bau','-append');
    
    %% BAU case
    disp('Calculating technology trajectory under business-as-usual case.');
    % BAU case
    % Set retirement to zero
    lbLDV = zeros(N_LDV*analysis_horizon*(Tcrit_LDV+1),1);
    ubLDV = Inf*ones(N_LDV*analysis_horizon*(Tcrit_LDV+1),1);
    
    lbLDVbau = lbLDV;
    ubLDVbau = ubLDV;
    % ubLDVbau(N_LDV*analysis_horizon+1:end) = 0;
    BdemLDVconstraint = BdemLDV + totfleetdemand;
    
    ALDV = [-AstockLDV;-Acafe;-AcleanLDV];
    BLDV = [-BstockLDV;-Bcafe;-BcleanLDV];
    AeqLDV = [AdemLDV];
    BeqLDV = [BdemLDVconstraint];
    
    options = optimset('Display','iter','LargeScale','on','MaxIter',1000,'TolFun',1e-3,'UseParallel','always');
    [xbau, fvalbau, exitflagbau, outputbau] = linprog(fLDV,ALDV,BLDV,AeqLDV,BeqLDV,lbLDVbau,ubLDVbau,[],options);
    
    Xbau(:,run) = xbau;
    % BAU trajectorry data
    Ebau = AemLDVbau*xbau - BemLDVbau;
    newbau = reshape(xbau(1:N_LDV*analysis_horizon),N_LDV,analysis_horizon);
    retbau = reshape(xbau(N_LDV*analysis_horizon+1:end),N_LDV,Tcrit_LDV,analysis_horizon);
    rettotbau = squeeze(sum(retbau,2));
    
    oldbau = AstockLDV*xbau-BstockLDV;
    oldbau = reshape(oldbau,N_LDV,Tcrit_LDV,analysis_horizon);
    oldtotbau = squeeze(sum(oldbau,2));
    
    totbau = oldtotbau + newbau;
    
    for k=1:analysis_horizon
        for i=1:N_LDV
            for j=2:Tcrit_LDV + 1
                if k==1
                    deadbau(i,j-1,k) = initfleetLDV(i,j-1).*disc_prob_LDV(i,j-1);
                else
                    deadbau(i,j-1,k) = oldbau(i,j-1,k-1).*disc_prob_LDV(i,j);
                end
            end
        end
    end
    
    deadtotbau = squeeze(sum(deadbau,2));
    
    for k=1:analysis_horizon
        mpgebau(k) = sum(FE_combined(:,k).*newbau(:,k),1)/sum(newbau(:,k),1);
    end
    
    % Total costs by year and type
    totLDVcostmatbau = [fLDVnewyearly*xbau, fLDVoldyearly*xbau - fLDVconstyearly, fLDVretyearly*xbau];
    
    % Electricity demand in MWh from EV charging
    EVchargedemandbau = zeros(analysis_horizon,1);
    for k=1:analysis_horizon
        for i=3:N_LDV
            EVchargedemandbau(k) = EVchargedemandbau(k) + Daily_Charging_Load(i,k)*totbau(i,k)*365/1000;
        end
    end
    EVCHARGEDEMBAU(:,run) = EVchargedemandbau;
    
    NEWTOTBAU(:,run) = sum(newbau(:,1:target_year-analysis_year0),2);
    OLDTOTBAU(:,run) = sum(oldtotbau(:,1:target_year-analysis_year0),2);
    RETTOTBAU(:,run) = sum(-rettotbau(:,1:target_year-analysis_year0),2);
    DEADTOTBAU(:,run) = sum(-deadtotbau(:,1:target_year-analysis_year0),2);
    EBAU(:,run) = Ebau(1:target_year-analysis_year0);
    
    % Write BAU data to spreadsheet
    xlswrite(excel_filename,newbau(:,1:target_year-analysis_year0)','BAU','B55');
    xlswrite(excel_filename,-deadtotbau(:,1:target_year-analysis_year0)','BAU','F55');
    xlswrite(excel_filename,-rettotbau(:,1:target_year-analysis_year0)','BAU','J55');
    xlswrite(excel_filename,totbau(:,1:target_year-analysis_year0)','BAU','N55');
    xlswrite(excel_filename,Ebau(1:target_year-analysis_year0),'BAU','S55');
    
    xlswrite(excel_filename,mpgebau(1:target_year-analysis_year0)','BAU','V55');
    xlswrite(excel_filename,FE_cafe(1:target_year-analysis_year0)','BAU','W55');
    
    xlswrite(excel_filename,EVchargedemandbau(1:target_year-analysis_year0),'BAU','Y55');
    
    %% Optimization under climate action case with delays
    % cay = climate action year
    disp('Calculating technology trajectory under climate action.');
    
    cuml_incrLDVcost = zeros(caylength,1);
    cuml_incrnewLDVcost = zeros(caylength,1);
    cuml_incroldLDVcost = zeros(caylength,1);
    cuml_incrretLDVcost = zeros(caylength,1);
    cuml_incrsocLDVcost = zeros(caylength,1);
    
    lbLDVmod = lbLDV;
    ubLDVmod = ubLDV;
    
    em_deficit = 0;
    linear_emcon_of_interest = linear_emcon((analysis_year0 - reference_year + 1) + 1:end);   %portion of ideal emission curve from 2011 to end of analysis time horizon
    
    for cay=1:caylength
        
        status = strcat('Climate Action Year: ',num2str(analysis_year0 + cay));
        disp(status);
        
        % Non-BAU case
        disp('Calculting technology trajectory under climate action case');
        
        %% Climate action emission factors and constraints setup
        disp('Setting up climate action emission factors and constraint.');
        elecCO2content = squeeze(elecCO2contentmat(cay,:));
        if cay>1
            elecCO2content(1:cay-1) = elecCO2contentbau(1:cay-1);  %Assign EFs from BAU to years before CAY
        end
        % Recompute A matrix and B vector for emission constraint based on
        % updated electric charging emission factors
        [AemLDV, BemLDV, EF_LDV_total_new, EF_LDV_total_old] =...
            EmissionConstraintSetupLDV(N_LDV, Tcrit_LDV, analysis_horizon, gasmiles, elecmiles, FE_gas, FE_elec, FEdecline_rate, AstockLDV, BstockLDV, elecCO2content);
        save(mat_filename,'AemLDV','BemLDV','EF_LDV_total_new','EF_LDV_total_old','-append');
        
        AemLDVconstraint(1,:) = sum(AemLDV(1:target_year - analysis_year0,:),1);
        BemLDVconstraint(1) = sum(BemLDV(1:target_year - analysis_year0)) + gross_emission_budget - emission_debt - cuml_emission_lockin;
        
        % Add a flat yearly emission constraint to years beyond target_year
        if analysis_horizon>target_year - analysis_year0
            AemLDVconstraint(2:analysis_horizon - (target_year - analysis_year0) + 1,:) = AemLDV(target_year - analysis_year0 + 1:analysis_horizon,:);
            BemLDVconstraint(2:analysis_horizon - (target_year - analysis_year0) + 1) = BemLDV(target_year - analysis_year0 + 1:analysis_horizon);
            
            for year=target_year + 1:analysis_year0 + analysis_horizon
                index = year - target_year + 1;
                BemLDVconstraint(index) = BemLDVconstraint(index) + (ideal_em_curve_slope*(target_year - reference_year) + emission_history(1));
            end
            BemLDVconstraint = BemLDVconstraint';
        end
        BemLDVconstraint_modified = BemLDVconstraint;   %initialize modified emission budget with the climate action in 2015 value
        
        % Modify emissions budget with any deficit
        disp('Setting up emission constraint.');
        BemLDVconstraint_modified(1) = BemLDVconstraint_modified(1) - em_deficit; % subtract accrued emission deficit from 2015 to climate action year from emissions budget
        
        %% Modify clean tech deployment constraint to exclude years before
        % current climate action year
        [AcleanLDVmod, BcleanLDVmod] = CleanTechConstraintSetup(N_LDV, Tcrit_LDV, analysis_horizon, 2015 - analysis_year0 + cay - 1, ref_index, dep_indices, ratios);
        
        %%
        if cay==4
            include_SCC = 1;
            SCC_inc_rate = 0.015;
            SCC_flag = 1;
            E_ref = 21.912e9;
            
            while SCC_flag > 0
                % Create new cost function with social cost of carbon incorporated
                disp('Adding carbon cost if applicable');
                fLDV_SCC = zeros(size(fLDV));
                fLDV_SCC_yearly = zeros(analysis_horizon,length(fLDV));
                fLDV_SCC_const_yearly = zeros(analysis_horizon,1);
                
                for k=1:analysis_horizon
                    fLDV_SCC = fLDV_SCC + (SCC*((1+SCC_inc_rate)^(k-1))/(1+discountrate)^(k-1))*AemLDV(k,:);    % 1.03 factor introduced because a 3% discount rate is built into the SCC values from literature
                    fLDV_SCC_yearly(k,:) = (SCC*((1+SCC_inc_rate)^(k-1))/(1+discountrate)^(k-1))*AemLDV(k,:);
                    fLDV_SCC_const_yearly(k) = (SCC*((1+SCC_inc_rate)^(k-1))/(1+discountrate)^(k-1))*BemLDV(k);
                end
                %     save('LDVMats.mat','fLDV_SCC','fLDV_SCC_yearly','fLDV_SCC_const_yearly','-append');
                
                %Prepare matrices for optimization
                
                ALDV = [-AstockLDV;-Acafe;-AcleanLDVmod];
                BLDV = [-BstockLDV;-Bcafe;-BcleanLDVmod];
                fLDV_effective = fLDV + fLDV_SCC;
                
                AeqLDV = [AdemLDV];
                BeqLDV = [BdemLDVconstraint];
                
                options = optimset('Display','iter','LargeScale','on','MaxIter',1000,'TolFun',1e-3,'UseParallel','always');
                
                [x, fval, exitflag, output] = linprog(fLDV_effective,ALDV,BLDV,AeqLDV,BeqLDV,lbLDVmod,ubLDVmod,[],options);
                
                E_SCC = AemLDV*x - BemLDV;
                E_check = sum(E_SCC(1:36));
                E_SCC_error = (E_check - E_ref)/E_ref;
                fprintf('\n\nEmissions error = %f percent\n',E_SCC_error*100);
                if abs(E_SCC_error) < 0.01
                    scc_prompt2 = '\n\nContinue? (Y/N = 1/0): ';
                    SCC_flag = input(scc_prompt2);
                    scc_prompt1 = 'Enter SCC value: ';
                    SCC = input(scc_prompt1);
                else
                    scc_prompt1 = 'Enter SCC value: ';
                    SCC = input(scc_prompt1);
                end
            end
        else
            ALDV = [-AstockLDV;-Acafe;-AcleanLDVmod;AemLDVconstraint];
            BLDV = [-BstockLDV;-Bcafe;-BcleanLDVmod;BemLDVconstraint_modified];
            fLDV_effective = fLDV;
            
            AeqLDV = [AdemLDV];
            BeqLDV = [BdemLDVconstraint];
            
            options = optimset('Display','iter','LargeScale','on','MaxIter',1000,'TolFun',1e-3,'UseParallel','always');
            
            [x, fval, exitflag, output] = linprog(fLDV_effective,ALDV,BLDV,AeqLDV,BeqLDV,lbLDVmod,ubLDVmod,[],options);
        end
        
        % Set emissions and technology trajectory in years prior to climate action to BAU
        % case for optimization in the following iteration
        for z=1:cay
            em_deficit = em_deficit + Ebau(z) - linear_emcon_of_interest(z);
        end
%         em_deficit = em_deficit + Ebau(cay) - linear_emcon_of_interest(cay);   % new emission deficit
        
        lbLDVmod(1:N_LDV*cay) = 0.95*xbau(1:N_LDV*cay);
        lbLDVmod(N_LDV*analysis_horizon + 1:N_LDV*analysis_horizon + N_LDV*cay*Tcrit_LDV) = 0.95*xbau(N_LDV*analysis_horizon + 1:N_LDV*analysis_horizon + N_LDV*cay*Tcrit_LDV);
        ubLDVmod(1:N_LDV*cay) = 1.05*xbau(1:N_LDV*cay);
        ubLDVmod(N_LDV*analysis_horizon + 1:N_LDV*analysis_horizon + N_LDV*cay*Tcrit_LDV) = 1.05*xbau(N_LDV*analysis_horizon + 1:N_LDV*analysis_horizon + N_LDV*cay*Tcrit_LDV);
        
        %% Calculate and store optimal fleet and emissions data
        disp('Calculating optimal fleet and emissions data.');
        
        new = reshape(x(1:N_LDV*analysis_horizon),N_LDV,analysis_horizon);
        ret = reshape(x(N_LDV*analysis_horizon+1:end),N_LDV,Tcrit_LDV,analysis_horizon);
        rettot = squeeze(sum(ret,2));
        
        old = AstockLDV*x-BstockLDV;
        old = reshape(old,N_LDV,Tcrit_LDV,analysis_horizon);
        oldtot = squeeze(sum(old,2));
        
        for k=1:analysis_horizon
            for i=1:N_LDV
                for j=2:Tcrit_LDV + 1
                    if k==1
                        dead(i,j-1,k) = initfleetLDV(i,j-1).*disc_prob_LDV(i,j-1);    % disc_prob age index is 'j' and not 'j-1' because the first colummn in it is for new vehicles (i.e. age = 0)
                    else
                        dead(i,j-1,k) = old(i,j-1,k-1).*disc_prob_LDV(i,j);
                    end
                end
            end
        end
        deadtot = squeeze(sum(dead,2));
        
        tot = oldtot + new;
        totyearly = sum(tot,1);
        for k=1:analysis_horizon
            totpercent(:,k) = tot(:,k)*100/totyearly(k);
        end
        
        retbyage = squeeze(sum(ret,1));
        for j=1:Tcrit_LDV
            for k=1:analysis_horizon
                temp(k,j) = retbyage(j,k);
            end
        end
        retbyage = temp;    % switch rows and columns for 3D bar graph such that color gradient displays change in age
        
        for k=1:analysis_horizon
            mpge(k) = sum(FE_combined(:,k).*new(:,k),1)/sum(new(:,k),1);
        end
        
        E = AemLDV*x - BemLDV;
        
        NEWTOT(:,cay,run) = sum(new(:,1:target_year-analysis_year0),2);
        OLDTOT(:,cay,run) = sum(oldtot(:,1:target_year-analysis_year0),2);
        RETTOT(:,cay,run) = sum(-rettot(:,1:target_year-analysis_year0),2);
        DEADTOT(:,cay,run) = sum(-deadtot(:,1:target_year-analysis_year0),2);
        ECA(:,cay,run) = E(1:target_year-analysis_year0);
        
        % Electricity demand in MWh from EV charging
        EVchargedemand = zeros(analysis_horizon,1);
        for k=1:analysis_horizon
            for i=3:N_LDV
                EVchargedemand(k) = EVchargedemand(k) + Daily_Charging_Load(i,k)*tot(i,k)*365/1000;
            end
        end
        EVCHARGEDEM(:,cay,run) = EVchargedemand;
        
        % Abatement costs by year and type
        incrnewLDVcosts = fLDVnewyearly*(x - xbau);
        incroldLDVcosts = fLDVoldyearly*(x - xbau);
        incrretLDVcosts = fLDVretyearly*(x - xbau);
        incrLDV_costs_by_type = [incrnewLDVcosts, incroldLDVcosts, incrretLDVcosts];
        incrLDVcostsyearly = sum(incrLDV_costs_by_type,2);
        
        % Total costs by year and type
        totLDVcostmat = [fLDVnewyearly*x, fLDVoldyearly*x - fLDVconstyearly, fLDVretyearly*x];
        
        % Ratio of total cost under climate action to total cost under bau
        Ctotaction = sum(sum(totLDVcostmat(1:target_year-analysis_year0,:)));
        Ctotbau = sum(sum(totLDVcostmatbau(1:target_year-analysis_year0,:)));
        R(cay,run) = Ctotaction/Ctotbau;
        CBAU(:,run) = Ctotbau;
        CACTION(cay,run) = Ctotaction;
        CINCR(cay,run) = sum(incrLDVcostsyearly(1:target_year-analysis_year0));
        CEFFECT(cay,run) = sum(incrLDVcostsyearly(1:target_year-analysis_year0))...
            /sum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0));
        
        X(:,cay,run) = x;
        cuml_incrLDVcost(cay) = sum(incrLDVcostsyearly(1:target_year-analysis_year0));
        cuml_incrnewLDVcost(cay) = sum(incrnewLDVcosts(1:target_year-analysis_year0));
        cuml_incroldLDVcost(cay) = sum(incroldLDVcosts(1:target_year-analysis_year0));
        cuml_incrretLDVcost(cay) = sum(incrretLDVcosts(1:target_year-analysis_year0));
        
        %% Write results to Excel files
        disp('Storing optimal fleet and emissions data.');
        tab = num2str(analysis_year0 + cay);
        if exitflag~=1
            xlswrite(excel_filename,{'Cannot achieve emission target!'},tab,'B1');
            cell = strcat('C',num2str(cay+1));
            R(cay:caylength,run) = 0;
            xlswrite(excel_filename,R(cay,run),'Summary',cell);
            break
        else
            xlswrite(excel_filename,sum(incrLDVcostsyearly(1:target_year-analysis_year0))/1e9,tab,'B1');
            xlswrite(excel_filename,sum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0))/1e9,tab,'B2');
            xlswrite(excel_filename,sum(incrLDVcostsyearly(1:target_year-analysis_year0))/sum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0)),tab,'B3');
            xlswrite(excel_filename,new(:,1:target_year-analysis_year0)',tab,'B55');
            xlswrite(excel_filename,-deadtot(:,1:target_year-analysis_year0)',tab,'F55');
            xlswrite(excel_filename,-rettot(:,1:target_year-analysis_year0)',tab,'J55');
            xlswrite(excel_filename,tot(:,1:target_year-analysis_year0)',tab,'N55');
            xlswrite(excel_filename,E(1:target_year-analysis_year0),tab,'R55');
            xlswrite(excel_filename,Ebau(1:target_year-analysis_year0),tab,'S55');
            xlswrite(excel_filename,linear_emcon_of_interest(1:target_year-analysis_year0),tab,'T55');
            
            xlswrite(excel_filename,mpge(1:target_year-analysis_year0)',tab,'V55');
            xlswrite(excel_filename,FE_cafe(1:target_year-analysis_year0)',tab,'W55');
            
            xlswrite(excel_filename,EVchargedemand(1:target_year-analysis_year0),tab,'Y55');
                    
            xlswrite(excel_filename,incrLDV_costs_by_type(1:target_year-analysis_year0,:),tab,'AA55');
            xlswrite(excel_filename,sum(incrLDV_costs_by_type(1:target_year-analysis_year0,:),2),tab,'AD55');
            xlswrite(excel_filename,cumsum(sum(incrLDV_costs_by_type(1:target_year-analysis_year0,:),2)),tab,'AE55');
            xlswrite(excel_filename,cumsum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0)),tab,'AF55');
            
            xlswrite(excel_filename,retbyage(1:target_year-analysis_year0,:),tab,'AH55');
            
            cell = strcat('C',num2str(cay+1));
            xlswrite(excel_filename,R(cay,run),'Summary',cell);
        end
        clear AemLDVconstraint BemLDVconstraint BemLDVconstraint_modified
    end
    clear AemLDVconstraint BemLDVconstraint BemLDVconstraint_modified
end
% save('.\MATS\LDV\Synthesis_LDV.mat','R','CBAU','CACTION','CINCR','CEFFECT','X','Xbau','NEWTOT','OLDTOT','RETTOT','DEADTOT','ECA',...
%     'NEWTOTBAU','OLDTOTBAU','RETTOTBAU','DEADTOTBAU','EBAU','EVCHARGEDEM');
toc
