close all
clear all

% Number of repetitions 
Nrep      = 100;
% If you want to use parallel pool
par       = true;
% Simulation time
Tmax      = 181;
% Reproductive basic number
R0        = 3;
% Incubation time in days
t_preinf  = 2;
% Contagious time in days
t_infect  = 6;
% Severe disease in days (binomial distribution)
t_disease = [3 1/4]; % R p
tdisease  = (1-t_disease(2))*t_disease(1)/t_disease(2);
frac(4) = 1 - nbincdf(t_preinf+t_infect+2,t_disease(1),t_disease(2));
frac(3) = 1 - nbincdf(t_preinf+t_infect/2+2,t_disease(1),t_disease(2));
frac(2) = 1 - nbincdf(t_preinf/2+2,t_disease(1),t_disease(2));
frac(1) = 1 - nbincdf(1,t_disease(1),t_disease(2));

% Number of individuals
Nind      = 1e5;
% Initial daily cases
n0 = 60;
% number of nursing homes per million
NH_rate   = 200;
% Extra death factor in nursing homes
NH_factor = 2;
% Probability that a worker infects inside the nursing home
P_treb    = 0.6;
% Probability that a resident infects inside the nurisng home
P_resi    = 0.8;
% Vaccination velocity as a % of daily population vaccinated
vac_vel   = 0.2;
% Restrictions profile can be automatic or fixed
restrictt = 'Automatic';
% Hospital capacity per million people
HLIM      = 1000;
% Minimum restrictions in %
Min_res   = 0;
Min_res = Min_res/100;
% Maximum restrictions in %
Max_res   = 80;
Max_res = Max_res/100;
% Country to use population pyramid
pop_struc = 'Europe';
pop       = xlsread('populations5.xlsx',pop_struc);
Prob_dead = xlsread('populations5.xlsx','Prob_dead');
Prob_hosp = xlsread('populations5.xlsx','Prob_hosp');
Prob_simp = xlsread('populations5.xlsx','Prob_Sympt');
pop_resi = xlsread('populations5.xlsx','Prob_NH');
% Time till first dose is effective
t1D_effec = 10;
% Time till second dose is effective
t2D_effec = 1;
% Porportion of nurses per resident
prop_nurs = 0.2;
p = sum(pop_resi(:).*pop(:));
P = sum(sum(pop(5:12,:)));
pop_treb = 0.*pop;
pop_treb(5:12,:) = p/P*prop_nurs;

% Vaccination strategy
strategy = 'Nursing homes + vulnerable';
% Posible vaccination strategies
%{'No vaccination','Nursing homes + age','Age','Vulnerable','Nursing homes + vulnerable','Contagious','Nursing homes + contagious','Random','Nursing homes + random'};

% Increase probability to live in a nursing home
fR = 1;
% Increase probability to die
fD = 1;
% Increase probability to be hospitalized
fH = 1;
pop_resi = pop_resi*fR;
Prob_dead = Prob_dead*fD;
Prob_hosp = Prob_hosp*fH;

pH = sum(Prob_hosp(:).*pop(:));
n0E = n0*2;
n0I = n0*6;
n0H = round(n0*9*pH);
tot_deads = 491316;
pD = sum(Prob_dead(:).*pop(:));
n0R = round(tot_deads/pD/747e6*Nind);
n0S = Nind - n0E - n0I - n0H - n0R;
N0 = [n0S n0E n0I n0H n0R];
HLIM = HLIM/1e6*Nind;
vac_vel = vac_vel/100*Nind;

% Vaccine configuration
% AstraZeneca
vaccine = 'AstraZeneca';
D2_effect = [0.30 0.70 0.90];
D1_effect = 0.95*D2_effect;

% Moderna
%vaccine = 'Moderna';
%D2_effect = [0.41 0.95 0.95];
%D1_effect = 0.84*D2_effect;

% Pfizer
%vaccine = 'Pfizer';
%D2_effect = [0.41 0.95 0.95];
%D1_effect = 0.55*D2_effect;

% Time between both doses
t_betw_do = 28;

list_var = {'D1_effect','D2_effect','HLIM','Max_res','Min_res','N0','NH_factor',...
    'NH_rate','Nind','Nrep','P_resi','P_treb','Prob_dead','Prob_hosp','Prob_simp',...
    'R0','Tmax','par','pop','pop_resi','pop_treb','prop_nurs','restrictt','strategy',...
    't1D_effec','t2D_effec','t_disease','t_infect','t_preinf','vac_vel','t_betw_do',...
    'vaccine','frac'};

clear parameters

noms = list_var;
for k = 1:length(noms)
    eval(['parameters.' noms{k} ' = ' noms{k} ';'])
end

[EVI,T] = IBM_Tmax(parameters);
save(['Data_00000.mat'],'EVI','T','parameters')

