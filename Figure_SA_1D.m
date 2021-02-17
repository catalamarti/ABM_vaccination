close all
clear all

R0_int   = [1.5 5];       %  1
R0_0     = 3;
HLIM_int = [25 100];      %  2
HLIM_0   = 50;
resM_int = [0.6 0.9];     %  3
resM_0   = 0.8;
n0_int   = [10 500];      %  4
n0_0     = 60;
t1De_int = [6.5 21.5];    %  5 int
t1De_0   = 10;
tpre_int = [1 6];         %  6
tpre_0   = 2;
tinf_int = [2 14];        %  7
tinf_0   = 6;
tdis_int = [5 20];        %  8
tdis_0   = 9;
PD_int   = [0.5 4];       %  9
PD_0     = 1;
PR_int   = [0.5 4];       % 10
PR_0     = 1;
PH_int   = [0.5 4];       % 11
PH_0     = 1;
Pres_int = [0.5 1];       % 12
Pres_0   = 0.8;
Ptre_int = [0.5 1];       % 13
Ptre_0   = 0.6;
vv_int   = [24.5 1000.5]; % 14 int
vv_0     = 200;
FH_int   = [0.5 2.7];     % 15
FH_0     = 2;
NH_int   = [100 400];     % 16
NH_0     = 200;
pN_int   = [0.1 0.3];     % 17
pN_0     = 0.2;
tbd_int  = [20.5 84.5];   % 18 int
tbd_0    = 21;
D1e_int  = [0.3 1];       % 19
D1e_0    = 0.78;
D2As_int = [0 0.95];      % 20
D2As_0   = 0.37;
D2Si_int = [0.55 1];      % 21
D2Si_0   = 0.87;
D2Se_int = [0.55 1];      % 22
D2Se_0   = 0.93;
PS_int   = [0.5 4];       % 23
PS_0     = 1;

Ns = 20;
variables = {'R0','HLIM','resM','n0','t1De','tpre','tinf','tdis','PD','PR',...
    'PH','Pres','Ptre','vv','FH','NH','pN','tbd','D1e','D2As','D2Si','D2Se',...
    'PS'};
vari      = {'R0','HLIM','Max_res','n0','t1D_effec','t_preinf','t_infect',...
    't_disease','fD','fR','fH','P_resi','P_treb','vac_vel','NH_factor',...
    'NH_rate','prop_nurs','t_betw_do','D1_effect(1)/parameters.D2_effect(1)',...
    'D2_effect(1)','D2_effect(2)','D2_effect(3)','fS'};
vari_lab  = {'R_0','H_{lim}','\lambda_{max}','n_0','t_{1D}','t_{inc}',...
    't_{con}','t_{sev}','f_D','f_R','f_H','P_{inf,nh,res}','P_{inf,nh,work}',...
    'r_{vac}','f_{NH}','r_{nh}','prop_{nh}','\Deltat_{dos}','D_{1,eff}','D_{2,asy}',...
    'D_{2,sym}','D_{2,sev}','f_S'};

nD = cell(length(variables),Ns);
nI = cell(length(variables),Ns);
nH = cell(length(variables),Ns);
nR = cell(length(variables),Ns);
nN = cell(length(variables),Ns);
nA = cell(length(variables),Ns);
xx = zeros(length(variables),Ns);

pD = 0.012179882378474;
pH = 0.026435699400882;
pR = 0.007284038562794;

for k = 1:23
    for i = 1:Ns
        load(['./Sensitivity_1D/Data_' num2str(k,'%02u') '_' num2str(i,'%02u') '.mat'])
        nD{k,i} = T(:,1)/1e3;
        nI{k,i} = T(:,4)*100;
        nH{k,i} = T(:,2)/1e3;
        nR{k,i} = T(:,13);
        nN{k,i} = T(:,12)./sum(parameters.pop(:).*parameters.pop_resi(:)*1e5)*100;
        nA{k,i} = T(:,10)*100;
        if strcmp(vari{k},'n0')
            xx(k,i) = parameters.N0(2)/2;
        elseif strcmp(vari{k},'t_disease')
            xx(k,i) = (1-parameters.t_disease(2))*parameters.t_disease(1)/parameters.t_disease(2);
        elseif k==9
            pd = sum(parameters.pop(:).*parameters.Prob_dead(:));
            xx(k,i) = pd/pD;
        elseif k==10
            pr = sum(parameters.pop(:).*parameters.pop_resi(:));
            xx(k,i) = pr/pR;
        elseif k==11
            ph = sum(parameters.pop(:).*parameters.Prob_hosp(:));
            xx(k,i) = ph/pH;
        elseif k==23
            eval(['lim_lim = ' variables{k} '_int;'])
            xx(k,i) = lim_lim(1) + (lim_lim(2)-lim_lim(1))*(i-1)/(Ns-1);
        else
            eval(['xx(k,i)=parameters.' vari{k} ';'])
        end
    end
end

C = lines(6);
va = {'nD','nH','nI','nN','nA','nR'};
ti = {'Cumulative deaths (%)','Cumulative hospitalized (%)',...
    'Cumulative infected (%)','Cumulative nh deaths (%)',...
    'Cumulative nh cases (%)','restrictions (days)'};

REF(1) = median(nD{5,6});
REF(2) = median(nH{5,6});
REF(3) = median(nI{5,6});
REF(4) = median(nN{5,6});
REF(5) = median(nA{5,6});
REF(6) = median(nR{5,6});

rx = 0.06;
mx = 0.06;
lx = 0.02;
dy = 0.1;
my = 0.1;
uy = 0.03;
DX = (1-2*mx-rx-lx)/3;
DY = (1-my-dy-uy)/2;
AX = [rx dy+my+DY; rx+mx+DX dy+my+DY; rx+2*mx+2*DX dy+my+DY;...
    rx dy; rx+mx+DX dy; rx+2*mx+2*DX dy];

for I = 1:23
    
    f = figure(1);
    clf
    f.Position = [0.0030    0.0463    1.1093    0.5927]*1e3;
    
    x = xx(I,:);
    
    for J = 1:6
        eval(['Y = ' va{J} ';'])
        ax = axes('units','normalized','Position',[AX(J,:) DX DY]);
        XX = [];
        YY = [];
        nn = 0;
        for k = 1:20
            if k==1 || xx(I,k-1)~=xx(I,k)
                XX = [XX;xx(I,k)*ones(length(Y{I,k}),1)];
                YY = [YY;Y{I,k}/REF(J)];
                nn = nn+1;
            end
        end
        plot([xx(I,1)-100 xx(I,20)+100],[1 1],'--k')
        hold on
        eval(['xy = ' variables{I} '_0;'])
        plot(xy*[1 1],[0 5],'--k')
        boxchart(XX,YY,'BoxWidth',0.1,'BoxFaceColor',C(J,:),'MarkerStyle','none')
        xxx = unique(xx(I,:));
        axis([2*xxx(1)-xxx(2) 2*xxx(end-1)-xxx(end) 0 2])
        XT = ax.XTick;
        ylabel(ti{J})
        xlabel(variables{I})
        
        cla
        w = 0.8/nn;
        xx_n = (xx(I,:)-xx(I,1))/(xx(I,end)-xx(I,1));
        xy_n = (xy-xx(I,1))/(xx(I,end)-xx(I,1));
        XX_n = (XX-xx(I,1))/(xx(I,end)-xx(I,1));
        XT_n = (XT-xx(I,1))/(xx(I,end)-xx(I,1));
        plot([-1 2],[1 1],'--k')
        hold on
        plot(xy_n*[1 1],[0 5],'--k')
        boxchart(XX_n,YY,'BoxWidth',w,'BoxFaceColor',C(J,:),'MarkerStyle','none')
        ax.XLim = [xx_n(1)-2*w xx_n(end)+2*w];
        wn = XT_n(2)-XT_n(1);
        xt = XT_n(1):wn:2;
        xt = [flip(XT_n(1):-wn:-1) xt(2:end)];
        xxt = xt*(xx(I,end)-xx(I,1))+xx(I,1);
        ax.XTick = xt(xxt>0.001);
        xtl = xxt(xxt>0.001);
        switch I
            case 3
                xtl = xtl*100;
            case 2
                xtl = xtl*10;
            case 14
                xtl = xtl/1000;
            case 19
                xtl = xtl*100;
            case 20
                xtl = xtl*100;
            case 21
                xtl = xtl*100;
            case 22
                xtl = xtl*100;
        end
        ax.XTickLabel = xtl;
        ax.YLim = [0 2];
        ylabel(ti{J})
        xlabel(vari_lab{I})
    end
    
    print(f,['./PNGs/App3_1D/Sensitivity_1D_' variables{I} '.png'],'-dpng','-r300')
    savefig(f,['./FIGs/App3_1D/Sensitivity_1D_' variables{I} '.fig'])
    
end