close all
clear all

vac = {'AstraZeneca','Pfizer','Moderna'};

f = dir('./Exploracio_D1e/*.mat');
N = length(f);

TBD = zeros(N,1);
VAC = zeros(N,1);
VSP = zeros(N,1);
D1E = zeros(N,1);

for k = 1:N
    load([f(k).folder '/' f(k).name],'parameters')
    VAC(k) = find(strcmp(vac,parameters.vaccine));
    TBD(k) = round(parameters.t_betw_do/7);
    VSP(k) = parameters.vac_vel;
    D1E(k) = round(parameters.D1_effect(1)/parameters.D2_effect(1)*100);
end

nam = {'_AZ_90','_Pf_55','_Mo_80','_Pf_90'};

S = [1 90; 2 50; 2 80; 2 90];

for ij = 1:length(nam)
    
    name_fig = ['/FIG_04' nam{ij}];
    vac = S(ij,1);
    tbd = 4:12;
    vsp = [25 50 100 200 400 1000];
    d1e = S(ij,2);
    REF = [1 3];
    
    Infected = zeros(length(tbd),length(vsp));
    Hospital = zeros(length(tbd),length(vsp));
    Deadssss = zeros(length(tbd),length(vsp));
    Restrict = zeros(length(tbd),length(vsp));
    
    for k = 1:length(tbd)
        for j = 1:length(vsp)
            i = find(VAC==vac & TBD==tbd(k) & VSP==vsp(j) & D1E==d1e);
            load([f(i).folder '/' f(i).name],'T')
            Infected(k,j) = mean(T(:,4));
            Hospital(k,j) = mean(T(:,2));
            Deadssss(k,j) = mean(T(:,1));
            Restrict(k,j) = mean(T(:,13));
        end
    end
    
    I1 = Infected/Infected(REF(1),REF(2))-1;
    I2 = Hospital/Hospital(REF(1),REF(2))-1;
    I3 = Deadssss/Deadssss(REF(1),REF(2))-1;
    I4 = Restrict/Restrict(REF(1),REF(2))-1;
    
    ff = figure;
    ff.Position = [25.6667   88.3333  730  546];
    lx = 0.07;
    mx = 0.05;
    rx = 0.03;
    DX = (1-lx-mx-rx)/2;
    ly = 0.1;
    my = 0.1;
    uy = 0.06;
    DY = (1-ly-my-uy)/2;
    
    lim1 = [-10 20];
    c1 = [0 1 0];
    c2 = [1 1 0];
    c3 = [1 0 0];
    N = 1e3;
    [COL1] = matrixcolors(c1,c2,c3,N,lim1);
    ax = axes('units','normalized','Position',[lx+0*(mx+DX) ly+1*(my+DY) DX DY]);
    imagesc(100*I1)
    for k = 1:size(I3,1)
        for j = 1:size(I3,2)
            if k==REF(1) && j==REF(2)
                text(j,k,'ref.','HorizontalAlignment','center','VerticalAlignment','middle')
            else
                text(j,k,num2str(round(100*I1(k,j)),'%+i'),'HorizontalAlignment','center','VerticalAlignment','middle')
            end
        end
    end
    ax.TickDir = 'out';
    ax.XTick = 1:6;
    ax.XTickLabel = vsp/1e5*200;
    ax.YTick = 1:9;
    ax.YTickLabel = tbd;
    colormap(ax,COL1);
    c=colorbar;
    caxis(lim1)
    %xlabel('Vaccination speed (% day^-^1)')
    ylabel('Dose interval (weeks)')
    title('Cumulative infections change (%)')
    text(-0.08,1.08,'A','Units','normalized','FontSize',16)
    
    lim = [-50 30];
    c1 = [0 1 0];
    c2 = [1 1 0];
    c3 = [1 0 0];
    N = 1e3;
    [COL] = matrixcolors(c1,c2,c3,N,lim);
    ax = axes('units','normalized','Position',[lx+1*(mx+DX) ly+1*(my+DY) DX DY]);
    imagesc(100*I2)
    for k = 1:size(I3,1)
        for j = 1:size(I3,2)
            if k==REF(1) && j==REF(2)
                text(j,k,'ref.','HorizontalAlignment','center','VerticalAlignment','middle')
            else
                text(j,k,num2str(round(100*I2(k,j)),'%+i'),'HorizontalAlignment','center','VerticalAlignment','middle')
            end
        end
    end
    ax.TickDir = 'out';
    ax.XTick = 1:6;
    ax.XTickLabel = vsp/1e5*200;
    ax.YTick = 1:9;
    ax.YTickLabel = tbd;
    colormap(ax,COL)
    c=colorbar;
    caxis(lim)
    %xlabel('Vaccination speed (% day^-^1)')
    %ylabel('Dose interval (weeks)')
    title('Cumulative hospitalizations change (%)')
    text(-0.08,1.08,'B','Units','normalized','FontSize',16)
    
    lim = [-35 50];
    c1 = [0 1 0];
    c2 = [1 1 0];
    c3 = [1 0 0];
    N = 1e3;
    [COL] = matrixcolors(c1,c2,c3,N,lim);
    ax = axes('units','normalized','Position',[lx+0*(mx+DX) ly+0*(my+DY) DX DY]);
    imagesc(100*I3)
    for k = 1:size(I3,1)
        for j = 1:size(I3,2)
            if k==REF(1) && j==REF(2)
                text(j,k,'ref.','HorizontalAlignment','center','VerticalAlignment','middle')
            else
                text(j,k,num2str(round(100*I3(k,j)),'%+i'),'HorizontalAlignment','center','VerticalAlignment','middle')
            end
        end
    end
    ax.TickDir = 'out';
    ax.XTick = 1:6;
    ax.XTickLabel = vsp/1e5*200;
    ax.YTick = 1:9;
    ax.YTickLabel = tbd;
    colormap(ax,COL)
    c=colorbar;
    caxis(lim)
    xlabel('Vaccination speed (% day^-^1)')
    ylabel('Dose interval (weeks)')
    title('Cumulative deaths change (%)')
    text(-0.08,1.08,'C','Units','normalized','FontSize',16)
    
    lim = [-75 10];
    c1 = [0 1 0];
    c2 = [1 1 0];
    c3 = [1 0 0];
    N = 1e3;
    [COL] = matrixcolors(c1,c2,c3,N,lim);
    ax = axes('units','normalized','Position',[lx+1*(mx+DX) ly+0*(my+DY) DX DY]);
    imagesc(100*I4)
    for k = 1:size(I3,1)
        for j = 1:size(I3,2)
            if k==REF(1) && j==REF(2)
                text(j,k,'ref.','HorizontalAlignment','center','VerticalAlignment','middle')
            else
                text(j,k,num2str(round(100*I4(k,j)),'%+i'),'HorizontalAlignment','center','VerticalAlignment','middle')
            end
        end
    end
    ax.TickDir = 'out';
    ax.XTick = 1:6;
    ax.XTickLabel = vsp/1e5*200;
    ax.YTick = 1:9;
    ax.YTickLabel = tbd;
    colormap(ax,COL)
    c=colorbar;
    caxis(lim)
    xlabel('Vaccination speed (% day^-^1)')
    %ylabel('Dose interval (weeks)')
    title('Cumulative restrictions change (%)')
    text(-0.08,1.08,'D','Units','normalized','FontSize',16)
    
    print(ff,['./PNGs' name_fig],'-dpng','-r600')
    savefig(ff,['./FIGs' name_fig])
    
end

function [COL] = matrixcolors(col1,col2,col3,N,lims)

x = linspace(lims(1),lims(2),N);

COL = zeros(N,3);
for k = 1:N
    if x(k)>0
        COL(k,:) = col2*(lims(2)-x(k))/lims(2)+col3*x(k)/lims(2);
    else
        COL(k,:) = col2*(lims(1)-x(k))/lims(1)+col1*x(k)/lims(1);
    end
end
COL = max(COL,0);
COL = min(COL,1);

end