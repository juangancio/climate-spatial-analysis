close all
clear all
clc


%% Plotting all times series


figure, set(gcf,'Position',[643 103 774 730])
tile=tiledlayout(2,1,'TileSpacing','compact');
tile.TileSpacing = 'compact';
tile.Padding = 'compact';

anom_sat=readmatrix('../final_ts/elnino_NOAA_monthly_anomaly.csv');
t_sat=datetime(1981,9:length(anom_sat)+8,1);

ONI=readmatrix('../ONI.csv');
t_oni=datetime(1950,1:length(ONI),1);

ONI=ONI(381:891);
t_oni=t_oni(381:891);

nexttile(1); hold on
yyaxis right
plot(t_oni(1:end),ONI(1:end),'k-','LineWidth',2)%,'Color',[.5,.5,.5])
set(gca,'FontSize',20,'YMinorTick','on')
grid on
grid minor
ylabel('Temp. Anom. (ºC)')


nexttile(2); hold on
yyaxis right
plot(t_oni(1:end),ONI(1:end),'k-','LineWidth',2)%,'Color',[.5,.5,.5])
set(gca,'FontSize',20,'YMinorTick','on')
grid on
grid minor
ylabel('Temp. Anom. (ºC)')



for lag=[1 2 4 8]
    nexttile(1);
yyaxis left
colororder('default')
MI=csvread(['../elnino_anom_MI_hor_L4_lag_' num2str(lag) '.csv']);

plot(t_sat,MI,'LineWidth',1), hold on
%plot(run_ave(MI,12),'LineWidth',1), hold on

set(gca,'FontSize',20,'YMinorTick','on')
grid on
grid minor
set(gca, 'SortMethod', 'depth')
nexttile(2);
yyaxis left
MI=csvread(['../elnino_anom_MI_ver_L4_lag_' num2str(lag) '.csv']);


%
plot(t_sat,MI,'LineWidth',1)
%plot(run_ave(MI,12),'LineWidth',1), 
hold on
set(gca,'FontSize',20,'YMinorTick','on')
grid on
grid minor

end

nexttile(1),%legend('spatial lag = 0.25º','spatial lag = 0.50º','spatial lag = 1.0º', ...
    %'spatial lag = 2.0º','sst anom.','location','northwest')
ylabel('hor. SMI (a.u.)')
set(gca,'YLim',[-.6,2])
yyaxis right
set(gca,'YLim',[-5,5])
fill([t_sat,fliplr(t_sat)],[zeros(size(anom_sat')),-5.*ones(size(anom_sat'))],'k','FaceAlpha',.2,'EdgeColor','none', 'HandleVisibility', 'off')

%yyaxis right
%plot(t_sat(1:end-2),anom_sat(1:end-2),'k-','LineWidth',1)
%annotation('arrow',[.5 .7],[.5 .7],'LineWidth',1.5)
%text(.5, .7,'lag','FontSize',20)

nexttile(2),legend('spatial lag = 0.25º','spatial lag = 0.50º','spatial lag = 1.0º', ...
   'spatial lag = 2.0º','sst anom.','location','northwest')
yyaxis left
ylabel('ver. SMI (a.u.)')
set(gca,'YLim',[-.6,2])
yyaxis right
set(gca,'YLim',[-5,5])
fill([t_sat,fliplr(t_sat)],[zeros(size(anom_sat')),-5.*ones(size(anom_sat'))],'k','FaceAlpha',.2,'EdgeColor','none','HandleVisibility','off')

% text(1,2,'a)','FontSize',20)
% text(1,2.2,'b)','FontSize',20)
% text(1,2.3,'c)','FontSize',20)
%% Significance testing

lag = 8;
L = 4;
figure, tiledlayout(3,1)

regs={'3','34','4'};
for i = 1:3
region=regs{i};

nexttile, hold on

MI=csvread(['mi_ts/elnino_anom_MI_hor_L' num2str(L) '_lag_' num2str(lag) '_region_' region '.csv']);

MI=MI(317:end);
nino_MI=MI(ONI(317:end)>0.5);
non_MI=MI(ONI(317:end)<=0.5);

%Test of unequal size and variance
[H,P] = ttest2(nino_MI,non_MI,'Vartype','unequal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MI=csvread(['../elnino_anom_MI_hor_L4_lag_' num2str(lag) '.csv']);
MI=csvread(['mi_ts/elnino_anom_MI_hor_L' num2str(L) '_lag_' num2str(lag) '_region_' region '.csv']);

ninos=[68,125,160,196,256,278,304,340,411,451,508];
m_ninos=[];
l_ninos=[];
j_ninos=[];
k_ninos=[];
m_pre=[];

m_ninos_cp=[];
m_ninos_ep=[];

m_pre_cp=[];
m_pre_ep=[];

co=[];
kinds={'EP','CP','CP','EP','CP','CP','EP','CP','EP','EP','EP?'};

for i =1:length(ninos)
    
    j=ninos(i);
    while ONI(j)>.5
        j=j-1;
    end
    
    j_co=ninos(i);
    while ONI(j_co)>0
        j_co=j_co-1;
    end


    if ninos(i)==508
        k=length(ONI);
    else
        k=ninos(i);
        while ONI(k)>.5
            k=k+1;
        end
    end
    l_ninos=[l_ninos k-j];
    k_ninos=[k_ninos k];
    j_ninos=[j_ninos j];
    co=[co j_co];
    m_ninos=[m_ninos mean(MI(j:k))];
    if i<11
    if kinds{i}=='EP'
        m_ninos_ep=[m_ninos_ep mean(MI(j:k))];
    elseif kinds{i}=='CP'
        m_ninos_cp=[m_ninos_cp mean(MI(j:k))];
    end
    end
end

%kinds={'EP','CP','CP','EP','CP','CP','EP','CP','EP','EP','EP?'};
years = {'1986', '1992','1994', '1997', '2002', '2004', '2006', '2009', '2015', ...
    '2018', '2023'};

% nexttile(1);
% for i =1:length(j_ninos)
% 
%         plot(t_oni(j_ninos(i):k_ninos(i)),ONI(j_ninos(i):k_ninos(i)),'r-','LineWidth',1.5)
%         text(t_oni(j_ninos(i)+5),4,kinds{i},'FontSize',18)
% 
% end

 %figure(10), hold on, grid on  ,% plot(t_oni,ONI,'k','LineWidth',2)
 %figure(11), hold on
   for i = 1:length(j_ninos)

        %figure(10),plot(t_oni(j_ninos(i):k_ninos(i)),ONI(j_ninos(i):k_ninos(i)),'r','LineWidth',2)
        
       %figure(10), plot(t_oni(2*j_ninos(i)-k_ninos(i):j_ninos(i)),ONI(2*j_ninos(i)-k_ninos(i):j_ninos(i)),'b','LineWidth',2)
        %same length, inmidiatly before
        
        %plot(t_oni(co(i)-(k_ninos(i)-j_ninos(i)):co(i)),ONI(co(i)-(k_ninos(i)-j_ninos(i)):co(i)),'b','LineWidth',1.5) 
        %same length, cool off period
   
        %plot(t_oni(co(i)-mean(l_ninos):co(i)),ONI(co(i)-mean(l_ninos):co(i)),'b','LineWidth',1.5) 
        %mean length, cool off period
        
        %plot(t_oni(j_ninos(i)-mean(l_ninos):j_ninos(i)),ONI(j_ninos(i)-mean(l_ninos):j_ninos(i)),'b','LineWidth',1.5)
        %mean length, inmidiatly after

        x1=MI(j_ninos(i):k_ninos(i)); %niños

        % x2=MI(2*j_ninos(i)-k_ninos(i):j_ninos(i));%same length, inmidiatly before
        % m_pre=[m_pre mean(MI(2*j_ninos(i)-k_ninos(i):j_ninos(i)))];

        x2=MI(j_ninos(i)-mean(l_ninos):j_ninos(i)); %mean length, inmidiatly before
        m_pre=[m_pre mean(MI(j_ninos(i)-mean(l_ninos):j_ninos(i)))];
        if i<11
            if kinds{i}=='EP'
                m_pre_ep=[m_pre_ep mean(MI(j_ninos(i)-mean(l_ninos):j_ninos(i)))];
            elseif kinds{i}=='CP'
                m_pre_cp=[m_pre_cp mean(MI(j_ninos(i)-mean(l_ninos):j_ninos(i)))];
            end
        end
        x = [x1; x2;];
        g = [ones(length(x1), 1); zeros(length(x2), 1);];
       %figure(11), 
       
       boxplot(x, g,'positions',[i-0.2,i+0.2],'widths',.3)
       text(i-0.2,0.05,kinds{i},'FontSize',18)
   end

    set(gca,'XTick',[1:11])
    set(gca,'XTickLabel',years)
    set(gca,'FontSize',18)
    set(gcf,'position',[455 420 897 334])
    xlabel('years')
    ylabel('hor. SMI (a.u.)')

% Test for mean MI on nino events vs non-nino anteceding periods

    [h,p]=ttest(m_ninos,m_pre);
    disp('p-value for mean MI during el nino:');
    disp(p);

    

% Test for mean MI on nino events vs non-nino anteceding periods (EP events)
    [h,p1]=ttest(m_ninos_ep,m_pre_ep);
    disp('p-value for mean MI during EP el nino:');
    disp(p1);

% Test for mean MI on nino events vs non-nino anteceding periods (CP events)
    [h,p2]=ttest(m_ninos_cp,m_pre_cp);
    disp('p-value for mean MI during CP el nino:');
    disp(p2);


        title(['region: ' region ', p-value (all) = ' num2str(p) ', p-value (EP) = ' num2str(p1) ', p-value (CP) = ' num2str(p2)])

   end
% figure(10),legend('ONI','El nino period selected','previous period selected','location','northoutside','orientation','horizontal')
% set(gca,'XTickLabel',years)
%     set(gca,'FontSize',18)
%     set(gcf,'position',[455 420 897 334])
%     ylabel('temperature anomaly')
%     xlabel('years')
% MI=csvread(['../elnino_anom_MI_ver_L4_lag_' num2str(lag) '.csv']);