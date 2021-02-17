function [EVI,T_m] = IBM_Tmax(parameters)

EVI = cell(parameters.Nrep,1);
T = cell(parameters.Nrep,1);

% Just to paralelize the code
if parameters.par
    parfor k = 1:parameters.Nrep
        [EVI{k},T{k}] = IBM_MC_residencies_0(parameters);
    end
else
    for k = 1:parameters.Nrep
        [EVI{k},T{k}] = IBM_MC_residencies_0(parameters);
    end
end

T_m = zeros(parameters.Nrep,13);
for k = 1:parameters.Nrep
    T_m(k,:) = T{k};
end

end
function [EVI,T] = IBM_MC_residencies_0(parameters)

names = fieldnames(parameters);
for k = 1:length(names)
    eval([names{k} '=parameters.' names{k} ';'])
end

npresi = Nind*NH_rate/1e6;

PS = Prob_simp;
PR = pop_resi;
PH = Prob_hosp./PS./(1-PR);
PT = pop_treb./(1-pop_resi);
PD = Prob_dead./PS./PH./(1-PR+PR*NH_factor);
fm = 1./PD;
fh = sqrt(NH_factor)*ones(size(PD));
fm = min(fh,fm);
fh = NH_factor./fm;
PHR = fh.*PH;
PDR = fm.*PD;

f = sum(Prob_hosp(:).*pop(:));
f1 = f*(1-D1_effect(3));
f2 = f*(1-D2_effect(3));

T = zeros(1,12);
EVI = zeros(Tmax,15);

nresi  = poissrnd(npresi);
while(nresi<npresi*0.5)
    nresi = poissrnd(npresi);
end

[conf,estat] = configuracio_inicial(parameters,PS,PH,PHR,PD,PDR,PR,PT,nresi);
IR = cell(1,nresi);
for k = 1:nresi
    IR{k} = find(conf(:,20) == k);
end

if strcmp(restrictt,'Automatic')
    res = zeros(Tmax,1);
    auto_res = true;
else
    res = restrictions_vect(fixed_res(:,2)',fixed_res(:,1)',Tmax);
    auto_res = false;
end

if auto_res
    res(1) = restriccions_auto(HLIM,estat,f,f1,f2,N0,Nind,Min_res,Max_res,frac,R0,t_infect);
end

for dia = 1:Tmax
    
    % Els susceptibles vacunats avancen
    estat(estat(:,1)==1 & conf(:,13)==dia,8) = 1;
    estat(estat(:,1)==1 & conf(:,25)==dia,8) = 2;
    
    % Els exposats avancen un dia
    estat(estat(:,2)==0 & estat(:,1)==2,1) = 3;
    estat(estat(:,1)==2,2) = estat(estat(:,1)==2,2) - 1;
    
    % Els infectats avancen un dia
    id = estat(:,3)==0 & estat(:,1)==3;
    estat(id & conf(:, 6)==1 & estat(:,8)==0,1) = 4;
    estat(id & conf(:, 6)==0 & estat(:,8)==0,1) = 5;
    estat(id & conf(:,17)==1 & estat(:,8)==1,1) = 4;
    estat(id & conf(:,17)==0 & estat(:,8)==1,1) = 5;
    estat(id & conf(:,23)==1 & estat(:,8)==2,1) = 4;
    estat(id & conf(:,23)==0 & estat(:,8)==2,1) = 5;
    
    estat(id & conf(:, 6)==1 & estat(:,8)==0,7) = 1;
    estat(id & conf(:,17)==1 & estat(:,8)==1,7) = 1;
    estat(id & conf(:,23)==1 & estat(:,8)==2,7) = 1;
    
    estat(estat(:,1)==3,3) = estat(estat(:,1)==3,3) - 1;
    
    % Els hospitalitzats avancen un dia
    id = estat(:,4)==0 & estat(:,1)==4;
    estat(id & conf(:, 7)==1 & estat(:,8)==0,1) = 6;
    estat(id & conf(:, 7)==0 & estat(:,8)==0,1) = 5;
    estat(id & conf(:,18)==1 & estat(:,8)==1,1) = 6;
    estat(id & conf(:,18)==0 & estat(:,8)==1,1) = 5;
    estat(id & conf(:,24)==1 & estat(:,8)==2,1) = 6;
    estat(id & conf(:,24)==0 & estat(:,8)==2,1) = 5;
    estat(estat(:,1)==4,4) = estat(estat(:,1)==4,4) - 1;
    
    % Noves infeccions
    n = 0;
    i_treb = find(estat(:,1)==3 & conf(:,19)==2);
    pinf = conf(i_treb,4)*(1-res(dia));
    n_treb = poissrnd(pinf);
    itreb = P_treb<rand(length(pinf),1);
    n = n + sum(n_treb(itreb));
    n_treb = n_treb(~itreb);
    i_treb = i_treb(~itreb);
    for k = 1:length(n_treb)
        if n_treb(k)>0
            id = randsample(IR{conf(i_treb(k),20)},n_treb(k),true);
            id = unique(id);
            estat(id(estat(id,1)==1 & estat(id,8)==0),1) = 2;
            estat(id(estat(id,1)==1 & estat(id,8)==1 & conf(id,15)==0),1) = 2;
            estat(id(estat(id,1)==1 & estat(id,8)==2 & conf(id,21)==0),1) = 2;
        end
    end
    i_resi = find(estat(:,1)==3 & conf(:,19)==1);
    pinf = conf(i_resi,4)*(1-res(dia));
    n_resi = poissrnd(pinf);
    iresi = P_resi<rand(length(pinf),1);
    n = n + sum(n_resi(iresi));
    n_resi = n_resi(~iresi);
    i_resi = i_resi(~iresi);
    for k = 1:length(n_resi)
        if n_resi(k)>0
            id = randsample(IR{conf(i_resi(k),20)},n_resi(k),true);
            id = unique(id);
            estat(id(estat(id,1)==1 & estat(id,8)==0),1) = 2;
            estat(id(estat(id,1)==1 & estat(id,8)==1 & conf(id,15)==0),1) = 2;
            estat(id(estat(id,1)==1 & estat(id,8)==2 & conf(id,21)==0),1) = 2;
        end
    end
    pinf = conf(estat(:,1)==3 & conf(:,19)==0,4)*(1-res(dia));
    n = n + sum(poissrnd(pinf));
    if n>0
        id = randsample(1:Nind,n,true);
        id = unique(id);
        estat(id(estat(id,1)==1 & estat(id,8)==0),1) = 2;
        estat(id(estat(id,1)==1 & estat(id,8)==1 & conf(id,15)==0),1) = 2;
        estat(id(estat(id,1)==1 & estat(id,8)==2 & conf(id,21)==0),1) = 2;
    end
    
    for k = [1:3 5:6]
        EVI(dia,k) = sum(estat(:,1)==k);
    end
    EVI(dia,4) = sum(estat(:,1)==4 & conf(:,19)~=1);
    EVI(dia,7) = (length(unique(conf(ismember(estat(:,1),[3 4]) & conf(:,20)>0,20))))/nresi;
    EVI(dia,8) = (length(unique(conf(ismember(estat(:,1),[3 4 5 6]) & conf(:,20)>0 & conf(:,11)==1,20))))/nresi;
    EVI(dia,9) = sum(ismember(estat(:,1),[3 4]) & conf(:,19)==2)/sum(conf(:,19)==2);
    EVI(dia,10) = sum(ismember(estat(:,1),[3 4 5 6]) & conf(:,19)==2 & conf(:,11)==1)/sum(conf(:,19)==2);
    EVI(dia,11) = sum(ismember(estat(:,1),[3 4]) & conf(:,19)==1)/sum(conf(:,19)==1);
    EVI(dia,12) = sum(ismember(estat(:,1),[3 4 5 6]) & conf(:,19)==1 & conf(:,11)==1)/sum(conf(:,19)==1);
    EVI(dia,14) = sum(estat(:,8)==1);
    EVI(dia,15) = sum(estat(:,8)==2);
    
    if auto_res
        res(dia+1) = restriccions_auto(HLIM,estat,f,f1,f2,EVI(dia,1:4),Nind,Min_res,Max_res,frac,R0,t_infect);
    end
    
end

conf(:,27) = estat(:,1);

[mh,th] = max(EVI(:,4));
T( 1) = sum(estat(:,1)==6);
T( 2) = sum(estat(:,7) & conf(:,19)~=1);
T( 3) = mh;
T( 4) = sum(estat(:,1)~=1 & conf(:,11)==1)/Nind;
T( 5) = th;
th = find(EVI(:,4)>0,1,'last')+1;
if isempty(th)
    th = nan;
end
T( 6) = th;
[~,ti] = max(EVI(:,3));
T( 7) = ti;
T( 8) = length(unique(conf(estat(:,1)~=1 & conf(:,19)>0 & conf(:,11)==1,20)))/nresi;
T( 9) = EVI(Tmax,10);
T(10) = EVI(Tmax,12);
[~,ta] = max(EVI(:,9)+EVI(:,11));
T(11) = ta;
T(12) = sum(estat(:,1)==6 & conf(:,19)==1);
T(13) = sum(res);

EVI(:,13) = res(1:size(EVI,1));

end
function [y] = generar_poissrnd(r,N)
y = poissrnd(r,[N 1]);
id = y==0;
while(sum(id)>0)
    y(id) = poissrnd(r,[sum(id) 1]);
    id = y==0;
end
end
function [y] = generar_nbinrnd(r,N)
y = nbinrnd(r(1),r(2),[N 1]);
id = y==0;
while(sum(id)>0)
    y(id) = nbinrnd(r(1),r(2),[sum(id) 1]);
    id = y==0;
end
end
function [conf,estat] = configuracio_inicial(parameters,PS,PH,PHR,PD,PDR,PR,PT,nresi)

names = fieldnames(parameters);
for k = 1:length(names)
    eval([names{k} '=parameters.' names{k} ';'])
end

efectivitat_vacuna1 = D1_effect;
efectivitat_vacuna2 = D2_effect;

conf = zeros(Nind,27);
pop = cumsum(pop(:));

r = rand(1,Nind);
age = sum(r>pop)+1;
id = age';
sex = 1*(age>19)+1;
age(age>19) = age(age>19)-19;
conf(:,1) = age';
conf(:,2) = sex';
conf(:,3) = 2*R0*rand(Nind,1);
conf(:,19) = PR(id)>rand(Nind,1);
ir = conf(:,19)==0;
conf(ir,19) = 2*(PT(id(ir))>rand(sum(ir),1));
conf(conf(:,19)==1,20) = randsample(1:nresi,sum(conf(:,19)==1),true);
pk = hist(conf(conf(:,20)>0,20),1:nresi);
conf(conf(:,19)==2,20) = randsample(1:nresi,sum(conf(:,19)==2),true,pk/sum(pk));

conf(:,5) = PS(id)>rand(Nind,1);

s = conf(:,5)==1 & conf(:,19)~=1;
conf(s,6) = PH(id(s))>rand(sum(s),1);
h = conf(:,6)==1 & conf(:,19)~=1;
conf(h,7) = PD(id(h))>rand(sum(h),1);

s = conf(:,5)==1 & conf(:,19)==1;
conf(s,6) = PHR(id(s))>rand(sum(s),1);
h = conf(:,6)==1 & conf(:,19)==1;
conf(h,7) = PDR(id(h))>rand(sum(h),1);

conf(:,8) = generar_poissrnd(t_preinf,Nind);
conf(:,9) = generar_poissrnd(t_infect,Nind);
conf(:,10) = generar_nbinrnd(t_disease,Nind);
conf(:,4) = conf(:,3)./conf(:,9);
conf(:,14) = PS(id).*PH(id);

Prob = zeros(size(PH,1)*size(PH,2),4);
Prob(:,2) = (1-PS(:))*(1-efectivitat_vacuna1(1));
Prob(:,3) = PS(:).*(1-PH(:)).*(1-efectivitat_vacuna1(2));
Prob(:,4) = PS(:).*PH(:)*(1-efectivitat_vacuna1(3));
Prob(:,1) = 1 - Prob(:,2) - Prob(:,3) - Prob(:,4);
Prob = cumsum(Prob,2);
P = Prob(id,:);
iv = sum(P < rand(Nind,1)*ones(1,4),2)+1;

conf(iv==1 & conf(:,19)~=1,15) = 1;
conf(iv>2  & conf(:,19)~=1,16) = 1;
conf(iv==4 & conf(:,19)~=1,17) = 1;
conf(iv==4 & conf(:,19)~=1 & PD(id)>rand(Nind,1),18) = 1;

Prob = zeros(size(PH,1)*size(PH,2),4);
Prob(:,2) = (1-PS(:))*(1-efectivitat_vacuna1(1));
Prob(:,3) = PS(:).*(1-PHR(:)).*(1-efectivitat_vacuna1(2));
Prob(:,4) = PS(:).*PHR(:)*(1-efectivitat_vacuna1(3));
Prob(:,1) = 1 - Prob(:,2) - Prob(:,3) - Prob(:,4);
Prob = cumsum(Prob,2);
P = Prob(id,:);
iv = sum(P < rand(Nind,1)*ones(1,4),2)+1;

conf(iv==1 & conf(:,19)==1,15) = 1;
conf(iv>2  & conf(:,19)==1,16) = 1;
conf(iv==4 & conf(:,19)==1,17) = 1;
conf(iv==4 & conf(:,19)==1 & PDR(id)>rand(Nind,1),18) = 1;

Prob = zeros(size(PH,1)*size(PH,2),4);
Prob(:,2) = (1-PS(:))*(1-efectivitat_vacuna2(1));
Prob(:,3) = PS(:).*(1-PH(:)).*(1-efectivitat_vacuna2(2));
Prob(:,4) = PS(:).*PH(:)*(1-efectivitat_vacuna2(3));
Prob(:,1) = 1 - Prob(:,2) - Prob(:,3) - Prob(:,4);
Prob = cumsum(Prob,2);
P = Prob(id,:);
iv = sum(P < rand(Nind,1)*ones(1,4),2)+1;

conf(iv==1 & conf(:,19)~=1,21) = 1;
conf(iv>2  & conf(:,19)~=1,22) = 1;
conf(iv==4 & conf(:,19)~=1,23) = 1;
conf(iv==4 & conf(:,19)~=1 & PD(id)>rand(Nind,1),24) = 1;

Prob = zeros(size(PH,1)*size(PH,2),4);
Prob(:,2) = (1-PS(:))*(1-efectivitat_vacuna2(1));
Prob(:,3) = PS(:).*(1-PHR(:)).*(1-efectivitat_vacuna2(2));
Prob(:,4) = PS(:).*PHR(:)*(1-efectivitat_vacuna2(3));
Prob(:,1) = 1 - Prob(:,2) - Prob(:,3) - Prob(:,4);
Prob = cumsum(Prob,2);
P = Prob(id,:);
iv = sum(P < rand(Nind,1)*ones(1,4),2)+1;

conf(iv==1 & conf(:,19)==1,21) = 1;
conf(iv>2  & conf(:,19)==1,22) = 1;
conf(iv==4 & conf(:,19)==1,23) = 1;
conf(iv==4 & conf(:,19)==1 & PDR(id)>rand(Nind,1),24) = 1;

if ~strcmp(strategy,'No vaccination')
    switch strategy
        case 'Random'
            risc = rand(Nind,1);
            NH_first = false;
        case 'Vulnerable'
            risc = conf(:,14);
            NH_first = false;
        case 'Contagious'
            risc = conf(:,3);
            NH_first = false;
        case 'Age'
            risc = conf(:,1);
            NH_first = false;
        case 'Nursing homes + random'
            risc = rand(Nind,1);
            NH_first = true;
        case 'Nursing homes + vulnerable'
            risc = conf(:,14);
            NH_first = true;
        case 'Nursing homes + contagious'
            risc = conf(:,3);
            NH_first = true;
        case 'Nursing homes + age'
            risc = conf(:,1);
            NH_first = true;
        otherwise
            error('no strategy')
    end
    
    if NH_first
        ir = find(conf(:,19)>0);
        [~,vac_order_r] = sort(risc(ir),'Descend');
        in = find(conf(:,19)==0);
        [~,vac_order_n] = sort(risc(in),'Descend');
        conf(:,26) = [ir(vac_order_r);in(vac_order_n)];
    else
        [~,vac_order] = sort(risc,'Descend');
        conf(:,26) = vac_order;
    end
    
    vac = [1 1];
    dia = 1;
    while (vac(1)<Nind || vac(2)<Nind)
        if dia<t_betw_do
            nvac = 2*vac_vel;
            id = vac(1)+(0:nvac-1);
            id = id(id<=Nind);
            ii = conf(id,26);
            conf(ii,13) = dia + t1D_effec;
            vac(1) = vac(1) + nvac;
        else
            if vac(1)<=Nind
                nvac = vac_vel;
                id = vac(1)+(0:nvac-1);
                id = id(id<=Nind);
                ii = conf(id,26);
                conf(ii,13) = dia + t1D_effec;
                vac(1) = vac(1) + nvac;
                id = vac(2)+(0:nvac-1);
                id = id(id<=Nind);
                ii = conf(id,26);
                conf(ii,25) = dia + t2D_effec;
                vac(2) = vac(2) + nvac;
            else
                nvac = 2*vac_vel;
                id = vac(2)+(0:nvac-1);
                id = id(id<=Nind);
                ii = conf(id,26);
                conf(ii,25) = dia + t2D_effec;
                vac(2) = vac(2) + nvac;
            end
        end
        dia = dia + 1;
    end
    
end

estat = zeros(Nind,8);
estat(:,2:4) = conf(:,8:10);
for k = 2:length(N0)
    if k == 4
        ii = randsample(find(estat(:,1)==0 & conf(:,6)==1 & conf(:,19)~=1),N0(k));
    else
        ii = randsample(find(estat(:,1)==0),N0(k));
    end
    estat(ii,1) = k;
    if k == 4
        estat(ii,7) = 1;
    end
    if k>1 && k<5 && ~isempty(ii)
        ii = ii(estat(ii,k)>1);
        n = unique(estat(ii,k));
        for i = 1:length(n)
            ij = ii(estat(ii,k)==n(i));
            estat(ij,k) = randi(n(i),length(ij),1);
        end
    end
end
estat(estat(:,1)==0) = 1;
conf(:,11) = estat(:,1);

end
function [res] = restrictions_vect(restriccions,dies_restriccions,Tmax)
tt = 1:Tmax;
res = zeros(Tmax,1);
res(dies_restriccions) = restriccions;
res(tt<=dies_restriccions(1)) = restriccions(1);
res(tt>=dies_restriccions(end)) = restriccions(end);
for i = 1:length(dies_restriccions)-1
    t = tt(tt>dies_restriccions(i) & tt<dies_restriccions(i+1));
    res(t) = restriccions(i) + (restriccions(i+1)-restriccions(i))*(t-dies_restriccions(i))/(dies_restriccions(i+1)-dies_restriccions(i));
end

end
function [res] = restriccions_auto(HLIM,estat,f,f1,f2,EVI,Nind,Min_res,Max_res,frac,R0,t_infect)

ff = f*sum(estat(:,8)==0)+f1*sum(estat(:,8)==1)+f2*sum(estat(:,8)==2);
ff = ff/Nind;
new = 2*R0/t_infect*EVI(1)/Nind*EVI(3)*ff;
PRED = EVI(4)*frac(4)+EVI(3)*frac(3)*ff+EVI(2)*frac(2)*ff;
res = 1 - (HLIM-PRED)/new;
res = min(res,Max_res);
res = max(res,Min_res);

end
function [x] = correct_parameters(pop_resi,Prob_dead,Prob_hosp,Prob_simp,fNH,fR,fD,fH)
PR = fR*pop_resi;
Ph = fH*Prob_hosp;
Pd = fH*Prob_dead;
PH = (fH*Prob_hosp)./(Prob_simp.*(1-fR*pop_resi));
PD = (fD*Prob_dead.*(1-fR*pop_resi))./(fH*Prob_hosp.*(1-fR*pop_resi+fNH*fR*pop_resi));
Pr = fNH*PD.*PH;

x = (sum(PD(:)>1) + sum(PD(:)<0) + sum(PH(:)>1) + sum(PH(:)<0) + ...
    sum(Ph(:)>1) + sum(Ph(:)<0) + sum(Pd(:)>1) + sum(Pd(:)<0) + ...
    sum(Pr(:)>1) + sum(Pr(:)<0) + sum(PR(:)>1) + sum(PR(:)<0)) == 0;

end