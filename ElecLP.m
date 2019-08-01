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

N_Util = 11;  % no. of vehicle technologies considered
Tcrit_Util = 60; % Critical age after which capacity is assumed to be scrapped
discountrate = 0.07;

reference_year = 2010;
target_year = 2050;
reduction_target = 0.71;    % %reduction in CO2 by target year
analysis_year0 = 2014;
analysis_yearend = target_year + Tcrit_Util;
analysis_horizon = analysis_yearend - analysis_year0;  % no. of years considered in LP analysis

include_SCC = 0;   % binary variable for whether or not to include social cost of carbon in the objective function; 1 = include, 0 = exclude
include_decom_costs = 0;    % binary variable for whether plant decommisioning costs should be included in retirement costs

%% Create design of experiments
disp('Creating design of experiments.');

caylength = 15;
param_num = 3;
exp_num = 3^param_num;
UtilExps = fullfact(3*ones(1,param_num))'-1;
CaseNames = {'Low','Nom','High'};
CaseNumbers = [0, 1, 2];
input_filename = 'Util_Master_Input_File.xlsx';

X = zeros(N_Util*analysis_horizon*(1 + Tcrit_Util),caylength,exp_num);
Xbau = zeros(N_Util*analysis_horizon*(1 + Tcrit_Util),exp_num);
R = zeros(caylength,exp_num);
S = zeros(caylength,exp_num);
CACTION = zeros(caylength,exp_num);
CBAU = zeros(1,exp_num);
CINCR = zeros(caylength,exp_num);
CEFFECT = zeros(caylength,exp_num);
NEWTOT = zeros(N_Util,caylength,exp_num);
OLDTOT = zeros(N_Util,caylength,exp_num);
RETTOT = zeros(N_Util,caylength,exp_num);
DEADTOT = zeros(N_Util,caylength,exp_num);
NEWTOTBAU = zeros(N_Util,exp_num);
OLDTOTBAU = zeros(N_Util,exp_num);
RETTOTBAU = zeros(N_Util,exp_num);
DEADTOTBAU = zeros(N_Util,exp_num);
ECA = zeros(target_year - analysis_year0,caylength,exp_num);
EBAU = zeros(target_year - analysis_year0,exp_num);
EVCHARGEEF = zeros(analysis_horizon,caylength,exp_num);
EVCHARGEEFBAU = zeros(analysis_horizon,exp_num);

%% Load common MAT files
load '.\MATS\Util\UtilMatsCommon.mat'

%% Setup generation stocks
disp('Loading fleet data from Excel files.');
% Initial (2014) fleet composition
initfleetUtil = xlsread(input_filename,'Nom','B3:BJ13'); % intial (2014) fleet of generation
inittotalUtil = sum(sum(initfleetUtil));
initnewUtil = initfleetUtil(:,1); % first row stores values for new generation in 2014
initoldUtil = initfleetUtil(:,2:end);   % old generation in 2014

% Survival probabilities
disc_prob_Util = xlsread(input_filename,'Nom','B17:BJ27'); % discard probabilities for Utils of different ages
P_Util = 1.-disc_prob_Util;   % survival probabilities for Utils of different ages
% 
% disp('Setting up generation capacity stocks.');
% [AstockUtil, BstockUtil] = StockSetupUtil(N_Util, Tcrit_Util, analysis_horizon, initfleetUtil, P_Util);
% save('.\MATS\Util\UtilMatsCommon.mat','AstockUtil','BstockUtil');

%% Annual demand constraint on total generation
% disp('Setting up demand constraint.');
% [AdemUtil, BdemUtil] = DemandConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, 'TotFleet');
% save('.\MATS\Util\UtilMatsCommon.mat','AdemUtil','BdemUtil','-append');

%% Emission constraint setup data
disp('Setting up data for emission constraint.');

CarbIntensity = xlsread(input_filename,'Parameters','D9:D19'); % will give a N_Util x 1 vector containing carbon intensity in kg CO2/mmBtu for all util technologies
Old_heatrates = xlsread(input_filename,'Nom','B157:BI167'); % will give a N_Util x T_crit matrix containing generation-weighted heat rates in Btu/kWh for initial fleet
% Compute A matrix and B vector for emission constraint
emission_history = xlsread('EmissionConstraints.xls','EmissionHistoryUtil');
em_2010 = emission_history(1); % CO2 emissions (tCO2) in 2010
reduction_percent_2050 = reduction_target;
cuml_emission_lockin = 0;   % this variable is to be used to account for cumulative emission lock-ins as estimated by the SD model

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

history_years = reference_year:analysis_year0;  % time array from 2005 to 2014
emission_debt = sum(emission_history - (ideal_em_curve_slope.*(history_years - reference_year) + emission_history(1))'); % Sum of emission deficits from reference year to initial year assuming reference case as climate action in reference year followed by a linear emission reduction till 2050

history_years = analysis_year0+1:target_year;
gross_emission_budget = sum(ideal_em_curve_slope.*(history_years - reference_year) + emission_history(1)); % sum allowable emission levels from analysis year 1 to target year
net_emission_budget = gross_emission_budget - emission_debt;

%% Market share constraint
% disp('Setting up market shares constraint.');
% % Markest shares of new capacity added
% MS = xlsread(input_filename,'Nom','B190:CS200');  % Projected market shares of new capacity additions
% ConstrainedTechbau = [1,3,6,7];  % Index of technologies to apply >= style market share constraints
% [AmsUtilbau, BmsUtilbau] = MSConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, MS, ConstrainedTechbau, 'TotFleet');
% save('.\MATS\Util\UtilMatsCommon.mat','AmsUtilbau','BmsUtilbau','-append');
%     ConstrainedTech = [5,6];  % Index of technologies to apply >= style market share constraints
%     [AmsUtil, BmsUtil] = MSConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, MS, ConstrainedTech, 'TotFleet');
%     save(mat_filename,'AmsUtil','BmsUtil','-append');

%% RPS constraint
% disp('Setting up RPS constraint.');
% RPS_value = 0.15;
% RPS_year = 2025;
% RenewableTech = [8,9,10,11];
% [ARPS, BRPS] = RPSSetupUtil(N_Util, Tcrit_Util, analysis_horizon, AstockUtil, BstockUtil, RPS_value, RenewableTech, RPS_year, analysis_year0);
% save('.\MATS\Util\UtilMatsCommon.mat','ARPS','BRPS','-append');

%% Run scenarios
for run=1:exp_num
    fprintf('\n\nRunning experiment # %d...\n\n',run)
    %% create output and mat filenames and write input data parameter values for experiment # run
    excel_filename = strcat('\Results\Util\CAD_Util_',num2str(run),'.xlsx');
    mat_filename = strcat('.\MATS\Util\CAD_Util_Mats_',num2str(run),'.mat');
    
    %% List of parameters for uncertainty analysis
    %
    % 1	Technology costs including fuel prices
    % 2	Heat rates and their decline with age
    % 3	Demand
    
    %% Load data that varies with scenario
%     load(mat_filename);
    disp('Loading cost data.');
    % Plant costs
    Util_techcosts = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(1,run))),'B115:CS125'); % will give a N_Util x Y matrix containaing levelized annual capital costs of each technology ($/MWh)
    Util_maintcosts = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(1,run))),'B129:CS139');    % will give a N_Util x Y matrix containaing total maintenance costs for tech ($/MWh)
    Util_decomcosts = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(1,run))),'B171:BI181');  % will give a N_Util x Tcrit_Util matrix of decommissioning costs
    service_life = xlsread(input_filename,'Parameters','C9:C19'); % will give a N_Util x 1 vector containing designed service life in years for all technologies
    loan_period = 20;   % capital financing period in years
    Util_fuelprices = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(1,run))),'B31:CS41');   % N_Util x Y matrix containing fuel price trends for all techs
    Util_CFs = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(1,run))),'B101:CS111');   % N_Util x Y matrix containing capacity factors for all techs
    Util_retcosts = zeros(N_Util,Tcrit_Util,analysis_horizon);
    Util_revenues = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(1,run))),'B143:CS153')*10; % will give a N_Util x Y matrix containing expected revenue in $/MWh
    
    for k=1:analysis_horizon
        for i=1:N_Util
            for j=1:Tcrit_Util
                if j>=k
                    Util_retcosts(i,j,k) = (((max(loan_period - j + 1,0)/loan_period)*Util_techcosts(i,1)) + (max(service_life(i) - j,0)*Util_revenues(i,k)/Util_CFs(i,k)))/(1+discountrate)^(k-1)...
                        + include_decom_costs*Util_decomcosts(i,j);  % first term is the loan liability, second term is the compensation for lost revenue, and third term is the decommisioning cost
                    % loan period costs not divided by discount rate since they
                    % are already in the form of annuity
                else
                    Util_retcosts(i,j,k) = (((max(loan_period - j + 1,0)/loan_period)*Util_techcosts(i,k-j)) + (max(service_life(i) - j,0)*Util_revenues(i,k)/Util_CFs(i,k)))/(1+discountrate)^(k-1)...
                        + include_decom_costs*Util_decomcosts(i,j);  % first term is the loan liability, second term is the compensation for lost revenue, and third term is the decommisioning cost
                end
            end
        end
    end
    % SCC = xlsread('Input_File_SA_Util.xlsx','SCC');   %Social cost of carbon in $/ton of CO2
    
    disp('Loading demand and market share data.');
    % Demand including EV charging
    demand = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(3,run))),'B185:CS185')';
    EVloadbau = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(3,run))),'B186:CS186')';
    EVloadcad = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(3,run))),'B219:CS233')';
    
    %% Emission constraint matrix setup
    disp('Loading emissions and heat rate data from Excel files.');
    
    % Technology emission factors
    New_heatrates = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(2,run))),'B87:CS97'); % will give a N_Util x Y matrix containing heat rates in Btu/kWh for new units by year
    FEdecline_rate = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(2,run))),'B216'); % Annual rate at which heat rates increase with every year of capacity ageing
    Util_heatrates = zeros(N_Util,Tcrit_Util+1,analysis_horizon);
    Util_heatrates(:,1,:) = New_heatrates;
    Util_heatrates(:,2:end,1) = Old_heatrates;
    
    disp('Setting up emission constraint matrices.');
    [AemUtil, BemUtil, EFnew_Util, EFold_Util] = EmissionConstraintSetupUtil(N_Util, Tcrit_Util, analysis_horizon,...
        Util_heatrates, FEdecline_rate, CarbIntensity, AstockUtil, BstockUtil);
    save(mat_filename,'AemUtil','BemUtil','EFnew_Util','EFold_Util');
    
    AemUtilconstraint(1,:) = sum(AemUtil(1:target_year - analysis_year0,:),1);
    BemUtilconstraint = zeros(analysis_horizon - (target_year - analysis_year0) + 1, 1);
    BemUtilconstraint(1) = sum(BemUtil(1:target_year - analysis_year0)) + gross_emission_budget - emission_debt - cuml_emission_lockin;
    
    % Add a flat yearly emission constraint to years beyond target_year
    if analysis_horizon>target_year - analysis_year0
        AemUtilconstraint(2:analysis_horizon - (target_year - analysis_year0) + 1,:) = AemUtil(target_year - analysis_year0 + 1:analysis_horizon,:);
        BemUtilconstraint(2:analysis_horizon - (target_year - analysis_year0) + 1) = BemUtil(target_year - analysis_year0 + 1:analysis_horizon);
        
        for year=target_year + 1:analysis_year0 + analysis_horizon
            index = year - target_year + 1;
            BemUtilconstraint(index) = BemUtilconstraint(index) + (ideal_em_curve_slope*(target_year - reference_year) + emission_history(1));
        end
    end
    

    %% Create constraint for proportionate renewables deployment
    disp('Creating renewables deployment constraint');
    ref_index = 8;  % wind is used as reference
    dep_indices = [9,10,11];    % solar PV, solar thermal, and geothermal expressed as fraction of wind deployment
    ratios = xlsread(input_filename,char(CaseNames(CaseNumbers==UtilExps(2,run))),'B212:B214');
    [AcleanUtil, BcleanUtil] = CleanTechConstraintSetup(N_Util, Tcrit_Util, analysis_horizon, 2015 - analysis_year0, ref_index, dep_indices, ratios);
    save(mat_filename,'AcleanUtil','BcleanUtil','-append');
    
    %% NPV objective function
    disp('Setting up objective function.');
    [fUtilnewyearly, fUtiloldyearly, fUtilretyearly, fUtilconstyearly, CFnew_Util, CFold_Util, CFret_Util] = NPVUtil(N_Util, Tcrit_Util, analysis_horizon, Util_heatrates, FEdecline_rate,...
        Util_techcosts, Util_maintcosts, Util_retcosts, Util_fuelprices, AstockUtil, BstockUtil);
    
    fUtil = zeros(1,N_Util*analysis_horizon*(1+Tcrit_Util));
    for k=1:analysis_horizon
        fUtilnewyearly(k,:) = fUtilnewyearly(k,:)/(1 + discountrate)^(k-1);
        fUtiloldyearly(k,:) = fUtiloldyearly(k,:)/(1 + discountrate)^(k-1);
        fUtilretyearly(k,:) = fUtilretyearly(k,:)/(1 + discountrate)^(k-1);
        fUtilconstyearly(k,:) = fUtilconstyearly(k,:)/(1 + discountrate)^(k-1);
        fUtil(1,:) = fUtil(1,:) + fUtilnewyearly(k,:) + fUtiloldyearly(k,:) + fUtilretyearly(k,:);
    end
    
    save(mat_filename, 'fUtilnewyearly', 'fUtiloldyearly', 'fUtilretyearly', 'fUtilconstyearly','fUtil',...
        'CFnew_Util', 'CFold_Util', 'CFret_Util','-append');
    
    % Create new cost function with social cost of carbon incorporated
%     disp('Creating cost functions under carbon price case.');
%     fUtil_SCC = zeros(size(fUtil));
%     fUtil_SCC_yearly = zeros(analysis_horizon,length(fUtil));
%     fUtil_SCC_const_yearly = zeros(analysis_horizon,1);
%     
%     SCC = 10;
%     for k=1:analysis_horizon
%         fUtil_SCC = fUtil_SCC + SCC*AemUtil(k,:);
%         fUtil_SCC_yearly(k,:) = SCC*AemUtil(k,:);
%         fUtil_SCC_const_yearly(k) = SCC*BemUtil(k);
%     end
%     save(mat_filename,'fUtil_SCC','fUtil_SCC_yearly','fUtil_SCC_const_yearly','-append');
    
    %% Optimization under BAU case
    disp('Calculating technology trajectory under business as usual case.');
    lbUtil = zeros(N_Util*analysis_horizon*(Tcrit_Util+1),1);
    ubUtil = Inf*ones(N_Util*analysis_horizon*(Tcrit_Util+1),1);
    
    % Set retirement to zero
    lbUtilbau = lbUtil;
    ubUtilbau = ubUtil;
%     ubUtilbau(N_Util*analysis_horizon+1:end) = 0;
    
    for k=1:analysis_horizon
%         ubUtilbau(squish(1,k,1,N_Util,Tcrit_Util)) = 0;    % set biomass to increase at no more than 1000 MW a year
        ubUtilbau(squish(6,k,1,N_Util,Tcrit_Util)) = 1000*365*24*0.898;   % limit new nuclear to no more than 1000 MW a year
        ubUtilbau(squish(7,k,1,N_Util,Tcrit_Util)) = 1000*365*24*0.42;   % limit new hydro to increase at no more than 1000 MW a year
    end
    
    BdemUtilconstraintbau = BdemUtil + demand + EVloadbau;
    linear_emcon_of_interest = linear_emcon((analysis_year0 - reference_year + 1) + 1:end);   %portion of ideal emission curve from 2011 to end of analysis time horizon
    
    AUtil = [-AstockUtil;-ARPS];
    BUtil = [-BstockUtil;-BRPS];
    AeqUtil = [AdemUtil;AcleanUtil];
    BeqUtil = [BdemUtilconstraintbau;BcleanUtil];
    
    options = optimset('Display','iter','LargeScale','on','MaxIter',1000,'TolFun',1e-3,'UseParallel','always');
    [xbau, fvalbau, exitflagbau, outputbau] = linprog(fUtil,AUtil,BUtil,AeqUtil,BeqUtil,lbUtilbau,ubUtilbau,[],options);
    
    Xbau(:,run) = xbau;
    % BAU trajectorry data
    Ebau = AemUtil*xbau - BemUtil;
    newbau = reshape(xbau(1:N_Util*analysis_horizon),N_Util,analysis_horizon);
    retbau = reshape(xbau(N_Util*analysis_horizon+1:end),N_Util,Tcrit_Util,analysis_horizon);
    rettotbau = squeeze(sum(retbau,2));
    retbyagebau = squeeze(sum(retbau,1));
    for j=1:Tcrit_Util
        for k=1:analysis_horizon
            temp1(k,j) = retbyagebau(j,k);
        end
    end
    retbyagebau = temp1;    % switch rows and columns for 3D bar graph such that color gradient displays change in age
    oldbau = AstockUtil*xbau - BstockUtil;
    oldbau = reshape(oldbau,N_Util,Tcrit_Util,analysis_horizon);
    oldtotbau = squeeze(sum(oldbau,2));
    totbau = oldtotbau + newbau;
    totyearlybau = sum(totbau,1);
    deadbau = zeros(N_Util, Tcrit_Util, analysis_horizon);
%     deadbau(:,1,:) = 0;
    for k=1:analysis_horizon
        for i=1:N_Util
            for j=2:Tcrit_Util + 1
                if k==1
                    deadbau(i,j-1,k) = initfleetUtil(i,j-1).*disc_prob_Util(i,j-1);
                else
                    deadbau(i,j-1,k) = oldbau(i,j-1,k-1).*disc_prob_Util(i,j);
                end
            end
        end
    end
    deadtotbau = squeeze(sum(deadbau,2));
    totUtilcostmatbau = [fUtilnewyearly*xbau, fUtiloldyearly*xbau - fUtilconstyearly, fUtilretyearly*xbau];
    
    % EV charging emission factors (t/MWh)
    EVcharegeEFsbau = Ebau./totyearlybau';
    EVCHARGEEFBAU(:,run) = EVcharegeEFsbau;
    
    NEWTOTBAU(:,run) = sum(newbau(:,1:target_year-analysis_year0),2);
    OLDTOTBAU(:,run) = sum(oldtotbau(:,1:target_year-analysis_year0),2);
    RETTOTBAU(:,run) = sum(-rettotbau(:,1:target_year-analysis_year0),2);
    DEADTOTBAU(:,run) = sum(-deadtotbau(:,1:target_year-analysis_year0),2);
    EBAU(:,run) = Ebau(1:target_year-analysis_year0);
    
    save(mat_filename,'NEWTOTBAU','OLDTOTBAU','RETTOTBAU','DEADTOTBAU','EBAU');
    
    xlswrite(excel_filename,newbau(:,1:target_year-analysis_year0)','BAU','B55');
    xlswrite(excel_filename,-deadtotbau(:,1:target_year-analysis_year0)','BAU','M55');
    xlswrite(excel_filename,-rettotbau(:,1:target_year-analysis_year0)','BAU','X55');
    xlswrite(excel_filename,totbau(:,1:target_year-analysis_year0)','BAU','AI55');
    xlswrite(excel_filename,Ebau(1:target_year-analysis_year0),'BAU','AU55');
    xlswrite(excel_filename,linear_emcon_of_interest(1:target_year-analysis_year0),'BAU','AV55');
    
    xlswrite(excel_filename,EVcharegeEFsbau(1:target_year-analysis_year0),'BAU','AX55');
    
    xlswrite(excel_filename,retbyagebau(1:target_year-analysis_year0,:),'BAU','BG55');
    
    %% Optimization under delayed climate action case
    % cay = climate action year
    disp('Calculating technology trajectory under climate action.');
    
    cuml_incrUtilcost = zeros(caylength,1);
    cuml_incrnewUtilcost = zeros(caylength,1);
    cuml_incroldUtilcost = zeros(caylength,1);
    cuml_incrretUtilcost = zeros(caylength,1);
    cuml_incrsocUtilcost = zeros(caylength,1);
    
    lbUtilmod = lbUtil;
    ubUtilmod = ubUtil;
%     ubUtilmod(N_Util*analysis_horizon+1:end) = 0;
    
    % Restrict generation of nuclear, hydro, geothermal, and biomass. Set CCS
    % to zero
    for k=1:analysis_horizon
        ubUtilmod(squish(1,k,1,N_Util,Tcrit_Util)) = 0;    % set coal to 0
%         ubUtilmod(squish(5,k,1,N_Util,Tcrit_Util)) = 1000*365*24;    % set biomass to increase at no more than 1000 MW a year;
        ubUtilmod(squish(6,k,1,N_Util,Tcrit_Util)) = 1000*365*24*0.898;   % limit new nuclear to no more than 1000 MW a year
        ubUtilmod(squish(7,k,1,N_Util,Tcrit_Util)) = 1000*365*24*0.42;   % limit new hydro to increase at no more than 1000 MW a year
    end
    
    em_deficit = 0;
    BemUtilconstraint_modified = BemUtilconstraint;   %initialize modified emission budget with the climate action in 2015 value
    
    for cay=1:caylength
        
        status = strcat('Climate Action Year: ',num2str(analysis_year0 + cay));
        disp(status);
        
        % Non-BAU case
        disp('Calculting technology trajectory under climate action case');
        
        % Evcharging demand
        BdemUtilconstraint = BdemUtil + demand + squeeze(EVloadcad(:,cay));
        
        % Emission constraint
        disp('Setting up emission constraint.');
        BemUtilconstraint_modified(1) = BemUtilconstraint_modified(1) - em_deficit; % subtract accrued emission deficit from 2011 to climate action year from emissions budget
        
        disp('Modifying renewables deployment constraint.');
        % Modify clean tech deployment constraint to exclude years before
        % current climate action year
        [AcleanUtilmod, BcleanUtilmod] = CleanTechConstraintSetup(N_Util, Tcrit_Util, analysis_horizon, 2015 - analysis_year0 + cay - 1, ref_index, dep_indices, ratios);
        
        % Prepare matrices for optimization
        if include_SCC==1
            AUtil = [-AstockUtil;-ARPS];
            BUtil = [-BstockUtil;-BRPS];
            fUtil_effective = fUtil + fUtil_SCC;
        else
            AUtil = [-AstockUtil;AemUtilconstraint;-ARPS];
            BUtil = [-BstockUtil;BemUtilconstraint_modified;-BRPS];
            fUtil_effective = fUtil;
        end
        
        AeqUtil = [AdemUtil;AcleanUtilmod];
        BeqUtil = [BdemUtilconstraint;BcleanUtilmod];
        
        options = optimset('Display','iter','LargeScale','on','MaxIter',1000,'TolFun',1e-3,'UseParallel','always');
        
        [x, fval, exitflag, output] = linprog(fUtil_effective,AUtil,BUtil,AeqUtil,BeqUtil,lbUtilmod,ubUtilmod,[],options);
        
        % Set emissions and technology trajectory in years prior to climate action to BAU
        % case for optimization in the following iteration
        em_deficit = em_deficit + Ebau(cay) - linear_emcon_of_interest(cay);   % new emission deficit
        
        lbUtilmod(1:N_Util*cay) = 0.9*xbau(1:N_Util*cay);
        lbUtilmod(N_Util*analysis_horizon + 1:N_Util*analysis_horizon + N_Util*cay*Tcrit_Util) = 0.9*xbau(N_Util*analysis_horizon + 1:N_Util*analysis_horizon + N_Util*cay*Tcrit_Util);
        ubUtilmod(1:N_Util*cay) = 1.1*xbau(1:N_Util*cay);
        ubUtilmod(N_Util*analysis_horizon + 1:N_Util*analysis_horizon + N_Util*cay*Tcrit_Util) = 1.1*xbau(N_Util*analysis_horizon + 1:N_Util*analysis_horizon + N_Util*cay*Tcrit_Util);
        
        %% Calculate and store optimal fleet and emissions data
        
        disp('Calculating and storing optimal fleet and emissions data.');
        
        new = reshape(x(1:N_Util*analysis_horizon),N_Util,analysis_horizon);
        ret = reshape(x(N_Util*analysis_horizon+1:end),N_Util,Tcrit_Util,analysis_horizon);
        rettot = squeeze(sum(ret,2));
        
        old = AstockUtil*x-BstockUtil;
        old = reshape(old,N_Util,Tcrit_Util,analysis_horizon);
        oldtot = squeeze(sum(old,2));
        
        dead = zeros(N_Util, Tcrit_Util, analysis_horizon);
        for k=1:analysis_horizon
            for i=1:N_Util
                for j=2:Tcrit_Util + 1
                    if k==1
                        dead(i,j-1,k) = initfleetUtil(i,j-1).*disc_prob_Util(i,j-1);
                    else
                        dead(i,j-1,k) = oldbau(i,j-1,k-1).*disc_prob_Util(i,j);
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
        for j=1:Tcrit_Util
            for k=1:analysis_horizon
                temp(k,j) = retbyage(j,k);
            end
        end
        retbyage = temp;    % switch rows and columns for 3D bar graph such that color gradient displays change in age
        
        E = AemUtil*x - BemUtil;
        
        NEWTOT(:,cay,run) = sum(new(:,1:target_year-analysis_year0),2);
        OLDTOT(:,cay,run) = sum(oldtot(:,1:target_year-analysis_year0),2);
        RETTOT(:,cay,run) = sum(-rettot(:,1:target_year-analysis_year0),2);
        DEADTOT(:,cay,run) = sum(-deadtot(:,1:target_year-analysis_year0),2);
        ECA(:,cay,run) = E(1:target_year-analysis_year0);
        
        % EV charging emission factors (t/MWh)
        EVchargeEFs = E./totyearly';
        EVCHARGEEF(:,cay,run) = EVchargeEFs;
        
        % Costs by year and type
        newUtilcosts = fUtilnewyearly*(x - xbau);
        oldUtilcosts = fUtiloldyearly*(x - xbau);
        retUtilcosts = fUtilretyearly*(x - xbau);
%         socUtilcosts = fUtil_SCC_yearly*(x - xbau);
        totUtilcostmat = [fUtilnewyearly*x, fUtiloldyearly*x - fUtilconstyearly, fUtilretyearly*x];   % emission social cost added with a '-' sign because emissions defined as E = Aem*x - Bem
        
        incrUtil_costs_by_type = [newUtilcosts, oldUtilcosts, retUtilcosts];
        incrUtilcostsyearly = sum(incrUtil_costs_by_type,2);
        
        % Ratio of total cost under climate action to total cost under bau
        Ctotaction = sum(sum(totUtilcostmat(1:target_year-analysis_year0,:)));
        Ctotbau = sum(sum(totUtilcostmatbau(1:target_year-analysis_year0,:)));
        R(cay,run) = Ctotaction/Ctotbau;
        CBAU(:,run) = Ctotbau;
        CACTION(cay,run) = Ctotaction;
        CINCR(cay,run) = sum(incrUtilcostsyearly(1:target_year-analysis_year0));
        CEFFECT(cay,run) = sum(incrUtilcostsyearly(1:target_year-analysis_year0))...
            /sum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0));
        
        X(:,cay,run) = x;
        cuml_incrUtilcost(cay) = sum(incrUtilcostsyearly(1:target_year-analysis_year0));
        cuml_incrnewUtilcost(cay) = sum(newUtilcosts(1:target_year-analysis_year0));
        cuml_incroldUtilcost(cay) = sum(oldUtilcosts(1:target_year-analysis_year0));
        cuml_incrretUtilcost(cay) = sum(retUtilcosts(1:target_year-analysis_year0));
%         cuml_incrsocUtilcost(cay) = sum(socUtilcosts(1:target_year-analysis_year0));
        
        %% Write results to Excel files
        disp('Writing results to spreadsheet.');
        tab = num2str(analysis_year0 + cay);
        if exitflag~=1
            xlswrite(excel_filename,{'Cannot achieve emission target!'},tab,'B1');
            cell = strcat('C',num2str(cay+1));
            R(cay:caylength,run) = 0;
            xlswrite(excel_filename,R(cay,run),'Summary',cell);
            break
        else
            xlswrite(excel_filename,sum(incrUtilcostsyearly(1:target_year-analysis_year0))/1e9,tab,'B1');
            xlswrite(excel_filename,sum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0))/1e9,tab,'B2');
            xlswrite(excel_filename,sum(incrUtilcostsyearly(1:target_year-analysis_year0))/sum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0)),tab,'B3');
            xlswrite(excel_filename,new(:,1:target_year-analysis_year0)',tab,'B55');
            xlswrite(excel_filename,-deadtot(:,1:target_year-analysis_year0)',tab,'M55');
            xlswrite(excel_filename,-rettot(:,1:target_year-analysis_year0)',tab,'X55');
            xlswrite(excel_filename,tot(:,1:target_year-analysis_year0)',tab,'AI55');
            xlswrite(excel_filename,E(1:target_year-analysis_year0),tab,'AT55');
            xlswrite(excel_filename,Ebau(1:target_year-analysis_year0),tab,'AU55');
            xlswrite(excel_filename,linear_emcon_of_interest(1:target_year-analysis_year0),tab,'AV55');
            
            xlswrite(excel_filename,EVchargeEFs(1:target_year-analysis_year0),tab,'AX55');
            
            xlswrite(excel_filename,incrUtil_costs_by_type(1:target_year-analysis_year0,:),tab,'AZ55');
            xlswrite(excel_filename,sum(incrUtil_costs_by_type(1:target_year-analysis_year0,:),2),tab,'BC55');
            xlswrite(excel_filename,cumsum(sum(incrUtilcostsyearly(1:target_year-analysis_year0),2)),tab,'BD55');
            xlswrite(excel_filename,cumsum(Ebau(1:target_year-analysis_year0)-E(1:target_year-analysis_year0)),tab,'BE55');
            
            xlswrite(excel_filename,retbyage(1:target_year-analysis_year0,:),tab,'BG55');
            
            
            cell = strcat('C',num2str(cay+1));
            xlswrite(excel_filename,R,'Summary',cell);
        end
    end
%     save(mat_filename,'NEWTOT','OLDTOT','RETTOT','DEADTOT','ECA','CBAU','CACTION','CINCR','CEFFECT');
end
save('.\MATS\Util\Synthesis_Util.mat','R','CBAU','CACTION','CINCR','CEFFECT','X','Xbau','NEWTOT','OLDTOT','RETTOT','DEADTOT','ECA',...
    'NEWTOTBAU','OLDTOTBAU','RETTOTBAU','DEADTOTBAU','EBAU','EVCHARGEEF','EVCHAGEEFBAU');
toc
