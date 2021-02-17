close all
clear all

stra = {'No vaccination','Nursing homes + age','Age','Vulnerable','Nursing homes + vulnerable','Contagious','Nursing homes + contagious','Random','Nursing homes + random'};
vac = {'AstraZeneca','Pfizer','Moderna'};

f = dir('./PF90/*.mat');
N = length(f);

STR = zeros(N,1);
TBD = zeros(N,1);
VAC = zeros(N,1);
HLI = zeros(N,1);
VSP = zeros(N,1);

for k = 1:N
    load([f(k).folder '/' f(k).name],'parameters')
    STR(k) = find(strcmp(stra,parameters.strategy));
    VAC(k) = find(strcmp(vac,parameters.vaccine));
    TBD(k) = parameters.t_betw_do>40;
    HLI(k) = parameters.HLIM;
    VSP(k) = parameters.vac_vel;
    create_evolucio([f(k).folder '/' f(k).name])
    create_evolucioM([f(k).folder '/' f(k).name])
end

stra = {'NoVac','NH+A','Age','Vuln','NH+V','Cont','NH+C','Rand','NH+R'};

nam = {'_Pf_90','_Pf_90_faster','_Pf_90_slower'};

S1 = [5 2 0 50 100; 5 2 0 50 200; 5 2 0 50 50];
S2 = [5 2 1 50 100; 5 2 1 50 200; 5 2 1 50 50];

for ij = 1:length(nam)
    
    name_fig = ['/FIG_03' nam{ij}];
    s1 = S1(ij,:); %strategy vac tbd hlim vacvel
    s2 = S2(ij,:); %strategy vac tbd hlim vacvel
    
    i1 = find(STR==s1(1) & VAC==s1(2) & TBD==s1(3) & HLI==s1(4) & VSP==s1(5));
    load([f(i1).folder '/' f(i1).name],'evolucioM','Tmax','EVI');
    ev1 = evolucioM;
    EV1 = EVI;
    time1 = 1:Tmax;
    i2 = find(STR==s2(1) & VAC==s2(2) & TBD==s2(3) & HLI==s2(4) & VSP==s2(5));
    load([f(i2).folder '/' f(i2).name],'evolucioM','Tmax','EVI');
    ev2 = evolucioM;
    EV2 = EVI;
    time2 = 1:Tmax;
    
    N = 100;
    
    s1 = 1;
    s2 = 2;
    stra = {'3w','12w'};
    
    Nind = 1e5;
    id = [1:6];
    for k = 1:size(EV1)
        EV1{k}(:,id) = EV1{k}(:,id)/Nind;
        EV2{k}(:,id) = EV2{k}(:,id)/Nind;
    end
    ev1(:,id,:) = ev1(:,id,:)/Nind;
    ev2(:,id,:) = ev2(:,id,:)/Nind;
    
    lw = 0.3;
    ltnv = '--';
    
    ff = figure;
    ff.Position = [25.6667   88.3333  871.3333  546];
    lx = 0.1;
    mx = 0.1;
    rx = 0.03;
    DX = (1-lx-mx-rx)/2;
    ly = 0.1;
    my = 0.1;
    uy = 0.03;
    DY = (1-ly-my-uy)/2;
    
    c1 = [255   0   0]/255;
    c3 = [126   0   0]/255;
    
    ax0 = axes('units','normalized','Position',[lx ly+my+DY DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,ltnv,'linewidth',1.5,'Color',c3)
    for k = 1:N
        plot(time1,100*EV1{k}(:,2)+100*EV1{k}(:,3)+100*EV1{k}(:,4),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,2)+100*EV2{k}(:,3)+100*EV2{k}(:,4),'linewidth',lw,'Color',0.5*c3+0.5)
    end
    plot(time1,100*ev1(:,2,1)+100*ev1(:,3,1)+100*ev1(:,4,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,2,1)+100*ev2(:,3,1)+100*ev2(:,4,1),ltnv,'linewidth',1.5,'Color',c3)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0 Tmax];
    ax.YLim(1) = 0;
    ylabel('Covid-19 prevalence (%)')
    xlabel('Time (days)')
    legend(stra{s1(1)},stra{s2(1)},'Location','best')
    text(-0.15,1.05,'A','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    c3 = [  0  96 192]/255*0.8;
    c1 = [  0 192 192]/255;
    
    ax0 = axes('units','normalized','Position',[lx+mx+DX ly+my+DY DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,ltnv,'linewidth',1.5,'Color',c3)
    for k = 1:N
        plot(time1,100*EV1{k}(:,4),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,4),'linewidth',lw,'Color',0.5*c3+0.5)
    end
    plot(time1,100*ev1(:,4,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,4,1),ltnv,'linewidth',1.5,'Color',c3)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0 Tmax];
    ax.YLim(1) = 0;
    ylabel('Hospitalized prevalence (%)')
    xlabel('Time (days)')
    legend(stra{s1(1)},stra{s2(1)},'Location','best')
    text(-0.15,1.05,'B','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    c2 = [ 0   255 160]/255;
    c4 = [  0  126 100]/255;
    
    ax0 = axes('units','normalized','Position',[lx ly DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c2)
    hold on
    plot(nan,nan,ltnv,'linewidth',1.5,'Color',c4)
    for k = 1:N
        plot(time1,100*EV1{k}(:,6),'linewidth',lw,'Color',0.5*c2+0.5)
        plot(time2,100*EV2{k}(:,6),'linewidth',lw,'Color',0.5*c4+0.5)
    end
    plot(time1,100*ev1(:,6,1),'linewidth',1.5,'Color',c2)
    plot(time2,100*ev2(:,6,1),ltnv,'linewidth',1.5,'Color',c4)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0 Tmax];
    ax.YLim(1) = 0;
    ylabel('Deaths incidence (%)')
    xlabel('Time (days)')
    legend(stra{s1(1)},stra{s2(1)},'Location','best')
    text(-0.15,1.05,'C','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    c1 = [  0  32  64]/255;
    c2 = [ 64   0   0]/255;
    
    ax0 = axes('units','normalized','Position',[lx+mx+DX ly DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,ltnv,'linewidth',1.5,'Color',c2)
    for k = 1:N
        plot(time1,100*EV1{k}(:,13),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,13),'linewidth',lw,'Color',0.5*c2+0.5)
    end
    plot(time1,100*ev1(:,13,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,13,1),ltnv,'linewidth',1.5,'Color',c2)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    axis([0 Tmax 0 100])
    ylabel('Restrictions (%)')
    xlabel('Time (days)')
    legend(stra{s1(1)},stra{s2(1)},'Location','best')
    text(-0.15,1.05,'D','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    print(ff,['./PNGs' name_fig],'-dpng','-r600')
    savefig(ff,['./FIGs' name_fig])
    
end

function [] = create_evolucio(filename)

load(filename)

N = length(EVI);

if ~exist('evolucio','var')
    
    m = 0;
    for k = 1:length(EVI)
        m = max(m,size(EVI{k},1));
    end
    Tmax = m;
    
    for k = 1:length(EVI)
        ID = sum(EVI{k}==0,2);
        id = find(ID==15,1,'first');
        if ~isempty(id)
            EVI{k}(id:m,:) = EVI{k}(id-1,:);
        else
            n = size(EVI{k},1);
            EVI{k}(n:m,:) = EVI{k}(n,:);
        end
    end
    
    evolucio = zeros(m,15,3);
    for k = 1:m
        for i = 1:15
            me = zeros(N,1);
            for j = 1:N
                me(j) = EVI{j}(k,i);
            end
            evolucio(k,i,1) = median(me);
            q = quantile(me,[0.05 0.95]);
            evolucio(k,i,2) = q(1);
            evolucio(k,i,3) = q(2);
        end
    end
    save(filename,'evolucio','Tmax','-append')
    
end

end
function [] = create_evolucioM(filename)

load(filename)

N = length(EVI);

if ~exist('evolucioM','var')
    
    m = Tmax;
    
    for k = 1:length(EVI)
        ID = sum(EVI{k}==0,2);
        id = find(ID==15,1,'first');
        if ~isempty(id)
            EVI{k}(id:m,:) = EVI{k}(id-1,:);
        else
            n = size(EVI{k},1);
            EVI{k}(n:m,:) = EVI{k}(n,:);
        end
    end
    
    evolucioM = zeros(m,15,3);
    for k = 1:m
        for i = 1:15
            me = zeros(N,1);
            for j = 1:N
                me(j) = EVI{j}(k,i);
            end
            evolucioM(k,i,1) = mean(me);
            q = std(me);
            evolucioM(k,i,2) = evolucioM(k,i,1) - q;
            evolucioM(k,i,3) = evolucioM(k,i,1) + q;
        end
    end
    save(filename,'evolucioM','-append')
    
end

end