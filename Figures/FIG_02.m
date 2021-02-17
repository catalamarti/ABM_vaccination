close all
clear all

stra = {'No vaccination','Nursing homes + age','Age','Vulnerable','Nursing homes + vulnerable','Contagious','Nursing homes + contagious','Random','Nursing homes + random'};
vac = {'AstraZeneca','Pfizer','Moderna'};

f = dir('./Different_Strategies/*.mat');
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

nam = {'','_faster','_slower','_higher','_smaller'};
S1 = [8 1 0 50 100; 8 1 0 50 200; 8 1 0 50 50; 8 1 0 100 100; 8 1 0 10 100];
S2 = [3 1 0 50 100; 3 1 0 50 200; 3 1 0 50 50; 3 1 0 100 100; 3 1 0 10 100];
S3 = [4 1 0 50 100; 4 1 0 50 200; 4 1 0 50 50; 4 1 0 100 100; 4 1 0 10 100];
S4 = [5 1 0 50 100; 5 1 0 50 200; 5 1 0 50 50; 5 1 0 100 100; 5 1 0 10 100];
S5 = [6 1 0 50 100; 6 1 0 50 200; 6 1 0 50 50; 6 1 0 100 100; 6 1 0 10 100];

lin = 'mean';

for ij = 1:length(nam)
    
    s1 = S1(ij,:); %strategy vac tbd hlim vacvel
    s2 = S2(ij,:); %strategy vac tbd hlim vacvel
    s3 = S3(ij,:); %strategy vac tbd hlim vacvel
    s4 = S4(ij,:); %strategy vac tbd hlim vacvel
    s5 = S5(ij,:); %strategy vac tbd hlim vacvel
    
    i1 = find(STR==s1(1) & VAC==s1(2) & TBD==s1(3) & HLI==s1(4) & VSP==s1(5));
    load([f(i1).folder '/' f(i1).name],'evolucio','Tmax','EVI','T','evolucioM');
    T1 = T;
    if strcmp(lin,'mean')
        ev1 = evolucioM;
    else
        ev1 = evolucio;
    end
    EV1 = EVI;
    time1 = 1:Tmax;
    i2 = find(STR==s2(1) & VAC==s2(2) & TBD==s2(3) & HLI==s2(4) & VSP==s2(5));
    load([f(i2).folder '/' f(i2).name],'evolucio','Tmax','EVI','T','evolucioM');
    T2 = T;
    if strcmp(lin,'mean')
        ev2 = evolucioM;
    else
        ev2 = evolucio;
    end
    EV2 = EVI;
    time2 = 1:Tmax;
    i3 = find(STR==s3(1) & VAC==s3(2) & TBD==s3(3) & HLI==s3(4) & VSP==s3(5));
    load([f(i3).folder '/' f(i3).name],'evolucio','Tmax','EVI','T','evolucioM');
    T3 = T;
    if strcmp(lin,'mean')
        ev3 = evolucioM;
    else
        ev3 = evolucio;
    end
    EV3 = EVI;
    time3 = 1:Tmax;
    i4 = find(STR==s4(1) & VAC==s4(2) & TBD==s4(3) & HLI==s4(4) & VSP==s4(5));
    load([f(i4).folder '/' f(i4).name],'evolucio','Tmax','EVI','T','evolucioM');
    if strcmp(lin,'mean')
        ev4 = evolucioM;
    else
        ev4 = evolucio;
    end
    T4 = T;
    EV4 = EVI;
    time4 = 1:Tmax;
    i5 = find(STR==s5(1) & VAC==s5(2) & TBD==s5(3) & HLI==s5(4) & VSP==s5(5));
    load([f(i5).folder '/' f(i5).name],'evolucio','Tmax','EVI','T','evolucioM');
    if strcmp(lin,'mean')
        ev5 = evolucioM;
    else
        ev5 = evolucio;
    end
    EV5 = EVI;
    time5 = 1:Tmax;
    T5 = T;
    
    N = 100;
    
    name_fig = ['/FIG_02' nam{ij}];
    stra = {'Rand','Age','Vul','NH+vul','Cont'};
    
    Nind = 1e5;
    id = [1:6];
    for k = 1:size(EV1)
        EV1{k}(:,id) = EV1{k}(:,id)/Nind;
        EV2{k}(:,id) = EV2{k}(:,id)/Nind;
        EV3{k}(:,id) = EV3{k}(:,id)/Nind;
        EV4{k}(:,id) = EV4{k}(:,id)/Nind;
        EV5{k}(:,id) = EV5{k}(:,id)/Nind;
    end
    ev1(:,id,:) = ev1(:,id,:)/Nind;
    ev2(:,id,:) = ev2(:,id,:)/Nind;
    ev3(:,id,:) = ev3(:,id,:)/Nind;
    ev4(:,id,:) = ev4(:,id,:)/Nind;
    ev5(:,id,:) = ev5(:,id,:)/Nind;
    
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
    
%     c1 = [0.91 0.17 0.03];
%     c2 = [0.92 0.18 0.80];
%     c3 = [0.13 0.14 0.96];
%     c4 = [0.01 0.97 0.02];
%     c5 = [0.99 0.98 0.09];

    c1 = [0.00 0.00 0.00];
    c2 = [0.84 0.12 0.84];
    c3 = [0.00 0.40 0.95];
    c4 = [0.20 0.00 0.70];
    c5 = [0.95 0.12 0.12];
    
    ax0 = axes('units','normalized','Position',[lx ly+my+DY DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,'linewidth',1.5,'Color',c2)
    plot(nan,nan,'linewidth',1.5,'Color',c3)
    plot(nan,nan,'linewidth',1.5,'Color',c4)
    plot(nan,nan,'linewidth',1.5,'Color',c5)
    for k = 1:N
        plot(time1,100*EV1{k}(:,2)+100*EV1{k}(:,3)+100*EV1{k}(:,4),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,2)+100*EV2{k}(:,3)+100*EV2{k}(:,4),'linewidth',lw,'Color',0.5*c2+0.5)
        plot(time3,100*EV3{k}(:,2)+100*EV3{k}(:,3)+100*EV3{k}(:,4),'linewidth',lw,'Color',0.5*c3+0.5)
        plot(time4,100*EV4{k}(:,2)+100*EV4{k}(:,3)+100*EV4{k}(:,4),'linewidth',lw,'Color',0.5*c4+0.5)
        plot(time5,100*EV5{k}(:,2)+100*EV5{k}(:,3)+100*EV5{k}(:,4),'linewidth',lw,'Color',0.5*c5+0.5)
    end
    plot(time1,100*ev1(:,2,1)+100*ev1(:,3,1)+100*ev1(:,4,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,2,1)+100*ev2(:,3,1)+100*ev2(:,4,1),'linewidth',1.5,'Color',c2)
    plot(time3,100*ev3(:,2,1)+100*ev3(:,3,1)+100*ev3(:,4,1),'linewidth',1.5,'Color',c3)
    plot(time4,100*ev4(:,2,1)+100*ev4(:,3,1)+100*ev4(:,4,1),'linewidth',1.5,'Color',c4)
    plot(time5,100*ev5(:,2,1)+100*ev5(:,3,1)+100*ev5(:,4,1),'linewidth',1.5,'Color',c5)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0 Tmax];
    ax.YLim(1) = 0;
    ylabel('Covid-19 prevalence (%)')
    xlabel('Time (days)')
    %legend(stra)
    text(-0.15,1.05,'A','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ax0 = axes('units','normalized','Position',[lx+mx+DX ly+my+DY DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,'linewidth',1.5,'Color',c2)
    plot(nan,nan,'linewidth',1.5,'Color',c3)
    plot(nan,nan,'linewidth',1.5,'Color',c4)
    plot(nan,nan,'linewidth',1.5,'Color',c5)
    for k = 1:N
        plot(time1,100*EV1{k}(:,4),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,4),'linewidth',lw,'Color',0.5*c2+0.5)
        plot(time3,100*EV3{k}(:,4),'linewidth',lw,'Color',0.5*c3+0.5)
        plot(time4,100*EV4{k}(:,4),'linewidth',lw,'Color',0.5*c4+0.5)
        plot(time5,100*EV5{k}(:,4),'linewidth',lw,'Color',0.5*c5+0.5)
    end
    plot(time1,100*ev1(:,4,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,4,1),'linewidth',1.5,'Color',c2)
    plot(time3,100*ev3(:,4,1),'linewidth',1.5,'Color',c3)
    plot(time4,100*ev4(:,4,1),'linewidth',1.5,'Color',c4)
    plot(time5,100*ev5(:,4,1),'linewidth',1.5,'Color',c5)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0 Tmax];
    ax.YLim(1) = 0;
    ylabel('Hospitalized prevalence (%)')
    xlabel('Time (days)')
    %legend(stra)
    text(-0.15,1.05,'B','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ax0 = axes('units','normalized','Position',[lx ly DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,'linewidth',1.5,'Color',c2)
    plot(nan,nan,'linewidth',1.5,'Color',c3)
    plot(nan,nan,'linewidth',1.5,'Color',c4)
    plot(nan,nan,'linewidth',1.5,'Color',c5)
    for k = 1:N
        plot(time1,100*EV1{k}(:,6),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,6),'linewidth',lw,'Color',0.5*c2+0.5)
        plot(time3,100*EV3{k}(:,6),'linewidth',lw,'Color',0.5*c3+0.5)
        plot(time4,100*EV4{k}(:,6),'linewidth',lw,'Color',0.5*c4+0.5)
        plot(time5,100*EV5{k}(:,6),'linewidth',lw,'Color',0.5*c5+0.5)
    end
    plot(time1,100*ev1(:,6,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,6,1),'linewidth',1.5,'Color',c2)
    plot(time3,100*ev3(:,6,1),'linewidth',1.5,'Color',c3)
    plot(time4,100*ev4(:,6,1),'linewidth',1.5,'Color',c4)
    plot(time5,100*ev5(:,6,1),'linewidth',1.5,'Color',c5)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0 Tmax];
    ax.YLim(1) = 0;
    ylabel('Deaths incidence (%)')
    xlabel('Time (days)')
    legend(stra,'Location','Northwest')
    text(-0.15,1.05,'C','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ax0 = axes('units','normalized','Position',[lx+mx+DX ly DX DY]);
    plot(nan,nan,'linewidth',1.5,'Color',c1)
    hold on
    plot(nan,nan,'linewidth',1.5,'Color',c2)
    plot(nan,nan,'linewidth',1.5,'Color',c3)
    plot(nan,nan,'linewidth',1.5,'Color',c4)
    plot(nan,nan,'linewidth',1.5,'Color',c5)
    for k = 1:N
        plot(time1,100*EV1{k}(:,13),'linewidth',lw,'Color',0.5*c1+0.5)
        plot(time2,100*EV2{k}(:,13),'linewidth',lw,'Color',0.5*c2+0.5)
        plot(time3,100*EV3{k}(:,13),'linewidth',lw,'Color',0.5*c3+0.5)
        plot(time4,100*EV4{k}(:,13),'linewidth',lw,'Color',0.5*c4+0.5)
        plot(time5,100*EV5{k}(:,13),'linewidth',lw,'Color',0.5*c5+0.5)
    end
    plot(time1,100*ev1(:,13,1),'linewidth',1.5,'Color',c1)
    plot(time2,100*ev2(:,13,1),'linewidth',1.5,'Color',c2)
    plot(time3,100*ev3(:,13,1),'linewidth',1.5,'Color',c3)
    plot(time4,100*ev4(:,13,1),'linewidth',1.5,'Color',c4)
    plot(time5,100*ev5(:,13,1),'linewidth',1.5,'Color',c5)
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    axis([0 Tmax 0 100])
    ylabel('Restrictions (%)')
    xlabel('Time (days)')
    %legend(stra)
    text(-0.15,1.05,'D','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ff2 = figure;
    ff2.Position = [25.6667   88.3333  871.3333  546];
    lx = 0.1;
    mx = 0.1;
    rx = 0.03;
    DX = (1-lx-mx-rx)/2;
    ly = 0.1;
    my = 0.1;
    uy = 0.03;
    DY = (1-ly-my-uy)/2;
    
    c1 = [0.00 0.00 0.00];
    c2 = [0.84 0.12 0.84];
    c3 = [0.00 0.40 0.95];
    c4 = [0.20 0.00 0.70];
    c5 = [0.95 0.12 0.12];
    
    ax0 = axes('units','normalized','Position',[lx ly+my+DY DX DY]);
    
    I = 4;
    med = 1:5;
    for k = 1:5
        eval(['T = T' num2str(k) ';'])
        med(k) = median(T(:,I));
    end    
    [s,i] = sort(med);
    for k = 1:5
        eval(['T = T' num2str(i(k)) '/med(1);'])
        eval(['c = c' num2str(i(k)) ';'])
        boxchart(k*ones(100,1),T(:,I),'BoxFaceColor',c,'MarkerStyle','none')
        hold on
    end
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0.5 5.5];
    ylabel('Covid-19 cumulative incidence')
    ax.XTick = 1:5;
    ax.XTickLabel = stra(i);
    %legend(stra)
    text(-0.19,1.05,'A','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ax0 = axes('units','normalized','Position',[lx+mx+DX ly+my+DY DX DY]);
    I = 2;
    med = 1:5;
    for k = 1:5
        eval(['T = T' num2str(k) ';'])
        med(k) = median(T(:,I));
    end    
    [s,i] = sort(med);
    for k = 1:5
        eval(['T = T' num2str(i(k)) '/med(1);'])
        eval(['c = c' num2str(i(k)) ';'])
        boxchart(k*ones(100,1),T(:,I),'BoxFaceColor',c,'MarkerStyle','none')
        hold on
    end
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0.5 5.5];
    ylabel('Covid-19 cumulative hospitalized')
    ax.XTick = 1:5;
    ax.XTickLabel = stra(i);
    text(-0.19,1.05,'B','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ax0 = axes('units','normalized','Position',[lx ly DX DY]);
    I = 1;
    med = 1:5;
    for k = 1:5
        eval(['T = T' num2str(k) ';'])
        med(k) = median(T(:,I));
    end    
    [s,i] = sort(med);
    for k = 1:5
        eval(['T = T' num2str(i(k)) '/med(1);'])
        eval(['c = c' num2str(i(k)) ';'])
        boxchart(k*ones(100,1),T(:,I),'BoxFaceColor',c,'MarkerStyle','none')
        hold on
    end
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0.5 5.5];
    ylabel('Covid-19 cumulative deads')
    ax.XTick = 1:5;
    ax.XTickLabel = stra(i);
    text(-0.19,1.05,'C','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    ax0 = axes('units','normalized','Position',[lx+mx+DX ly DX DY]);
    I = 13;
    med = 1:5;
    for k = 1:5
        eval(['T = T' num2str(k) ';'])
        med(k) = median(T(:,I));
    end    
    [s,i] = sort(med);
    for k = 1:5
        eval(['T = T' num2str(i(k)) '/med(1);'])
        eval(['c = c' num2str(i(k)) ';'])
        boxchart(k*ones(100,1),T(:,I),'BoxFaceColor',c,'MarkerStyle','none')
        hold on
    end
    
    ax=gca;
    ax.LineWidth = 1.2;
    ax.Box = 'off';
    ax.TickDir = 'out';
    ax.XLim = [0.5 5.5];
    ylabel('Restrictions')
    ax.XTick = 1:5;
    ax.XTickLabel = stra(i);
    text(-0.19,1.05,'D','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
    
    print(ff,['./PNGs' name_fig],'-dpng','-r600')
    savefig(ff,['./FIGs' name_fig])
    print(ff2,['./PNGs' name_fig '_new'],'-dpng','-r600')
    savefig(ff2,['./FIGs' name_fig '_new'])
    
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