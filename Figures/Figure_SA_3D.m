close all
clear all

load('./Sensitivity_all/LHS.mat')

Np = size(LHS,2);
Ns = 1e3;

% EVI_tot = zeros(181,13,Ns);
% T_tot = zeros(Ns,13);
% 
% for k = 1:Ns
%     create_evolucio(['./Sensitivity_all/Data_' num2str(k,'%04u') '.mat'])
%     load(['./Sensitivity_all/Data_' num2str(k,'%04u') '.mat'],'T','evolucio','Tmax')
%     EVI_tot(:,:,k) = evolucio(:,1:13,1);
%     T_tot(k,:) = mean(T);
% end
% 
% save('./Sensitivity_all/LHS.mat','EVI_tot','T_tot','-append')
% 
% PRCC_EVI = zeros(181,13,Np);
% PRCC_T = zeros(13,Np);
% pPRCC_EVI = zeros(181,13,Np);
% pPRCC_T = zeros(13,Np);
% for k = 1:Np
%     X = LHS(:,k);
%     Z = LHS; Z(:,k)=[];
%     for i = 1:181
%         for j = 1:13
%             [PRCC_EVI(i,j,k),pPRCC_EVI(i,j,k)] = partialcorr(X,squeeze(EVI_tot(i,j,:)),Z);
%         end
%     end
%     for i = 1:13
%         [PRCC_T(i,k),pPRCC_T(i,k)] = partialcorr(X,T_tot(:,i),Z);
%     end
%     k
% end
% save('./Sensitivity_all/LHS.mat','PRCC_EVI','PRCC_T','-append')

tit = {'Total deaths','Total hospitalized','Maximum hospitalized','Total infected',...
    'Day peak hospital','last restrictions day','Peak of infected',...
    'NH affectation','workers affectation','residents affectation','NH peak',...
    'deaths in NH','Retsrictions'};
params = {'R_0','H_{lim}','\lambda_{max}','n_0','t_{1D}','t_{inc}',...
    't_{con}','t_{sev}','f_D','f_R','f_H','P_{inf,nh,res}','P_{inf,nh,work}',...
    'r_{vac}','f_{NH}','r_{nh}','prop_{nh}','\Deltat_{dos}','D_{1,eff}','D_{2,asy}',...
    'D_{2,sym}','D_{2,sev}','f_S'};

for k = 1:13
    [s,i] = sort(PRCC_T(k,:));
    figure
    bar(1:23,s)
    hold on
    plot([0 24],+0.082*[1 1],'--k')
    plot([0 24],-0.082*[1 1],'--k')
    axis([0.5 23.5 -1 1])
    ax=gca;
    ax.XTick = 1:23;
    ax.XTickLabel = params(i);
    ax.TickDir = 'out';
    xtickangle(90)
    ylabel('Partial correlation')
    title(tit{k})    
end

tit2 = {'Susceptible','Exposed','Infected','Hospitalized','Recovered','Deaths'};
tit2{13} = 'Restrictions';

for k = [1:6 13]
    figure
    P = squeeze(PRCC_EVI(:,k,:));
    ik = find(max(abs(P))>0.2);
    c = jet(length(ik));
    for i = 1:length(ik)
        plot(P(:,ik(i)),'linewidth',1.3,'Color',c(i,:))
        hold on
    end
    plot([0 240],+0.082*[1 1],'--k')
    plot([0 240],-0.082*[1 1],'--k')
    legend(params(ik),'Location','eastoutside')
    axis([0 180 -1 1])
    ylabel('Partial correlation')
    title(tit2{k})
end

for k = 1:23
    figure
    P = squeeze(PRCC_EVI(:,[1:6 13],k));
    c = jet(7);
    for i = 1:7
        plot(P(:,i),'linewidth',1.3,'Color',c(i,:))
        hold on
    end
    plot([0 240],+0.082*[1 1],'--k')
    plot([0 240],-0.082*[1 1],'--k')
    legend([tit2(1:6) tit2(13)],'Location','eastoutside')
    axis([0 180 -1 1])
    ylabel('Partial correlation')
    title(params{k})
end

params = {'R0','HLIM','resM','n0','t1De','tpre','tinf','tdis','PD','PR',...
    'PH','Pres','Ptre','vv','FH','NH','pN','tbd','D1e','D2As','D2Si','D2Se',...
    'PS'};

tit = [tit tit2(1:6) tit2(end) params];
for k = 1:43
    f = figure(k);
    print(f,['./PNGs/App3_all/Figure_' num2str(k,'%02u') '_' tit{k} '.png'],'-dpng','-r600')
    savefig(f,['./FIGs/App3_all/Figure_' num2str(k,'%02u') '_' tit{k} '.fig'])
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