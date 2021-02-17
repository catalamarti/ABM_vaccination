close all
clear all

stra = {'No vaccination','Nursing homes + age','Age','Vulnerable','Nursing homes + vulnerable','Contagious','Nursing homes + contagious','Random','Nursing homes + random'};
vac = {'AstraZeneca','Pfizer','Moderna'};

f = dir('./Different_countries/*.mat');
N = length(f);

countries = {f.name};
for k = 1:length(countries)
    i = find(countries{k}=='_');
    countries{k} = countries{k}(i(1)+1:i(end)-1);
end
count = unique(countries);

STR = zeros(N,1);
TBD = zeros(N,1);
VAC = zeros(N,1);
CNT = zeros(N,1);

for k = 1:N
    load([f(k).folder '/' f(k).name],'parameters')
    STR(k) = find(strcmp(stra,parameters.strategy));
    VAC(k) = find(strcmp(vac,parameters.vaccine));
    TBD(k) = parameters.t_betw_do>40;
    CNT(k) = find(strcmp(count,countries{k}));
    create_evolucio([f(k).folder '/' f(k).name])
end

stra = {'NoVac','NH+A','Age','Vuln','NH+V','Cont','NH+C','Rand','NH+R'};

for ijk = 1:length(count)
    
    nam = {'_AZ','_Pf','_Mo'};
    S1 = [1 1 0; 1 2 0; 1 3 0];
    S2 = [5 1 0; 5 2 0; 5 3 0];
    S3 = [5 1 1; 5 2 1; 5 3 1];
    
    for ij = 1:length(nam)
        
        s1 = S1(ij,:); %strategy vac tbd hlim vacvel
        s2 = S2(ij,:); %strategy vac tbd hlim vacvel
        s3 = S3(ij,:); %strategy vac tbd hlim vacvel
        
        i1 = find(STR==s1(1) & VAC==s1(2) & TBD==s1(3) & CNT==ijk);
        load([f(i1).folder '/' f(i1).name],'evolucio','Tmax','EVI');
        ev1 = evolucio;
        EV1 = EVI;
        time1 = 1:Tmax;
        i2 = find(STR==s2(1) & VAC==s2(2) & TBD==s2(3) & CNT==ijk);
        load([f(i2).folder '/' f(i2).name],'evolucio','Tmax','EVI');
        ev2 = evolucio;
        EV2 = EVI;
        time2 = 1:Tmax;
        i3 = find(STR==s3(1) & VAC==s3(2) & TBD==s3(3) & CNT==ijk);
        load([f(i3).folder '/' f(i3).name],'evolucio','Tmax','EVI');
        ev3 = evolucio;
        EV3 = EVI;
        time3 = 1:Tmax;
        
        N = 100;
        
        name_fig = ['/FIG_01_' count{ijk} nam{ij}];
        s1 = 1;
        s2 = 2;
        s3 = 3;
        stra = {'No vac.','NH+vul 3w','NH+vul 12w'};
        
        Nind = 1e5;
        id = [1:6];
        for k = 1:size(EV1)
            EV1{k}(:,id) = EV1{k}(:,id)/Nind;
            EV2{k}(:,id) = EV2{k}(:,id)/Nind;
            EV3{k}(:,id) = EV3{k}(:,id)/Nind;
        end
        ev1(:,id,:) = ev1(:,id,:)/Nind;
        ev2(:,id,:) = ev2(:,id,:)/Nind;
        ev3(:,id,:) = ev3(:,id,:)/Nind;
        
        lw = 0.3;
        ltnv = '--';
        ltvn = '-.';
        
        ff = figure(1);
        clf
        
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
        c2 = [126   0   0]/255;
        c3 = [126 126   0]/255;
        
        ax0 = axes('units','normalized','Position',[lx ly+my+DY DX DY]);
        plot(nan,nan,'linewidth',1.5,'Color',c1)
        hold on
        plot(nan,nan,ltnv,'linewidth',1.5,'Color',c2)
        plot(nan,nan,ltvn,'linewidth',1.5,'Color',c3)
        for k = 1:N
            plot(time1,100*EV1{k}(:,2)+100*EV1{k}(:,3)+100*EV1{k}(:,4),'linewidth',lw,'Color',0.5*c1+0.5)
            plot(time2,100*EV2{k}(:,2)+100*EV2{k}(:,3)+100*EV2{k}(:,4),'linewidth',lw,'Color',0.5*c2+0.5)
            plot(time3,100*EV3{k}(:,2)+100*EV3{k}(:,3)+100*EV3{k}(:,4),'linewidth',lw,'Color',0.5*c3+0.5)
        end
        plot(time1,100*ev1(:,2,1)+100*ev1(:,3,1)+100*ev1(:,4,1),'linewidth',1.5,'Color',c1)
        plot(time2,100*ev2(:,2,1)+100*ev2(:,3,1)+100*ev2(:,4,1),ltnv,'linewidth',1.5,'Color',c2)
        plot(time3,100*ev3(:,2,1)+100*ev3(:,3,1)+100*ev3(:,4,1),ltvn,'linewidth',1.5,'Color',c3)
        
        ax=gca;
        ax.LineWidth = 1.2;
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.XLim = [0 Tmax];
        ax.YLim(1) = 0;
        ylabel('Covid-19 prevalence (%)')
        xlabel('Time (days)')
        legend(stra,'Location','best')
        text(-0.15,1.05,'A','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
        
        c2 = [  0  96 192]/255*0.8;
        c1 = [  0 192 192]/255;
        c3 = [  0 192  96]/255;
        
        ax0 = axes('units','normalized','Position',[lx+mx+DX ly+my+DY DX DY]);
        plot(nan,nan,'linewidth',1.5,'Color',c1)
        hold on
        plot(nan,nan,ltnv,'linewidth',1.5,'Color',c2)
        plot(nan,nan,ltvn,'linewidth',1.5,'Color',c3)
        for k = 1:N
            plot(time1,100*EV1{k}(:,4),'linewidth',lw,'Color',0.5*c1+0.5)
            plot(time2,100*EV2{k}(:,4),'linewidth',lw,'Color',0.5*c2+0.5)
            plot(time3,100*EV3{k}(:,4),'linewidth',lw,'Color',0.5*c3+0.5)
        end
        plot(time1,100*ev1(:,4,1),'linewidth',1.5,'Color',c1)
        plot(time2,100*ev2(:,4,1),ltnv,'linewidth',1.5,'Color',c2)
        plot(time3,100*ev3(:,4,1),ltvn,'linewidth',1.5,'Color',c3)
        
        ax=gca;
        ax.LineWidth = 1.2;
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.XLim = [0 Tmax];
        ax.YLim(1) = 0;
        ylabel('Hospitalized prevalence (%)')
        xlabel('Time (days)')
        legend(stra,'Location','best')
        text(-0.15,1.05,'B','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
        
        c1 = [ 0   255 160]/255;
        c2 = [  0  126 100]/255;
        c3 = [  0   60 255]/255;
        
        ax0 = axes('units','normalized','Position',[lx ly DX DY]);
        plot(nan,nan,'linewidth',1.5,'Color',c1)
        hold on
        plot(nan,nan,ltnv,'linewidth',1.5,'Color',c2)
        plot(nan,nan,ltvn,'linewidth',1.5,'Color',c3)
        for k = 1:N
            plot(time1,100*EV1{k}(:,6),'linewidth',lw,'Color',0.5*c1+0.5)
            plot(time2,100*EV2{k}(:,6),'linewidth',lw,'Color',0.5*c2+0.5)
            plot(time3,100*EV3{k}(:,6),'linewidth',lw,'Color',0.5*c3+0.5)
        end
        plot(time1,100*ev1(:,6,1),'linewidth',1.5,'Color',c1)
        plot(time2,100*ev2(:,6,1),ltnv,'linewidth',1.5,'Color',c2)
        plot(time3,100*ev3(:,6,1),ltvn,'linewidth',1.5,'Color',c3)
        
        ax=gca;
        ax.LineWidth = 1.2;
        ax.Box = 'off';
        ax.TickDir = 'out';
        ax.XLim = [0 Tmax];
        ax.YLim(1) = 0;
        ylabel('Deaths incidence (%)')
        xlabel('Time (days)')
        legend(stra,'Location','best')
        text(-0.15,1.05,'C','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
        
        c1 = [  0   0  64]/255;
        c2 = [ 64   0   0]/255;
        c3 = [  0  64   0]/255;
        
        ax0 = axes('units','normalized','Position',[lx+mx+DX ly DX DY]);
        plot(nan,nan,'linewidth',1.5,'Color',c1)
        hold on
        plot(nan,nan,ltnv,'linewidth',1.5,'Color',c2)
        plot(nan,nan,ltvn,'linewidth',1.5,'Color',c3)
        for k = 1:N
            plot(time1,100*EV1{k}(:,13),'linewidth',lw,'Color',0.5*c1+0.5)
            plot(time2,100*EV2{k}(:,13),'linewidth',lw,'Color',0.5*c2+0.5)
            plot(time3,100*EV3{k}(:,13),'linewidth',lw,'Color',0.5*c3+0.5)
        end
        plot(time1,100*ev1(:,13,1),'linewidth',1.5,'Color',c1)
        plot(time2,100*ev2(:,13,1),ltnv,'linewidth',1.5,'Color',c2)
        plot(time3,100*ev3(:,13,1),ltvn,'linewidth',1.5,'Color',c3)
        
        ax=gca;
        ax.LineWidth = 1.2;
        ax.Box = 'off';
        ax.TickDir = 'out';
        axis([0 Tmax 0 100])
        ylabel('Restrictions (%)')
        xlabel('Time (days)')
        legend(stra,'Location','best')
        text(-0.15,1.05,'D','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontSize',16)
        
        print(ff,['./PNGs/App2' name_fig],'-dpng','-r600')
        savefig(ff,['./FIGs/App2' name_fig])
        
    end
    
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