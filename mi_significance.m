close all
clear all
clc


%% Plotting all times series



figure, set(gcf,'Position',[643 103 774 730])
tile=tiledlayout(3,1,'TileSpacing','compact');
tile.TileSpacing = 'compact';
tile.Padding = 'compact';

anom_sat=readmatrix('final_ts/elnino_NOAA_monthly_anomaly.csv');
t_sat=datetime(1981,9:length(anom_sat)+8,1);

ONI=readmatrix('input_data/ONI.csv');
t_oni=datetime(1950,1:length(ONI),1);

ONI=ONI(381:891);
t_oni=t_oni(381:891);

% nexttile(1); hold on
% yyaxis right
% plot(t_oni(1:end),ONI(1:end),'k-','LineWidth',2)%,'Color',[.5,.5,.5])
% set(gca,'FontSize',20,'YMinorTick','on')
% grid on
% grid minor
% ylabel('Temp. Anom. (ºC)')
% 
% 
% nexttile(2); hold on
% yyaxis right
% plot(t_oni(1:end),ONI(1:end),'k-','LineWidth',2)%,'Color',[.5,.5,.5])
% set(gca,'FontSize',20,'YMinorTick','on')
% grid on
% grid minor
% ylabel('Temp. Anom. (ºC)')

L = 4;
region = 'gulf';

for lag=[1 2 4 8]
    nexttile(1);
    %yyaxis left
    colororder('default')
    MI=csvread(['mi_ts/elnino_anom_MI_hor_L' num2str(L) '_lag_' num2str(lag) '_region_' region '.csv']);
    
    plot(t_sat,MI,'LineWidth',1.5), hold on
    %plot(run_ave(MI,12),'LineWidth',1), hold on
    
    set(gca,'FontSize',16,'YMinorTick','on')
    grid on
    grid minor
    set(gca, 'SortMethod', 'depth')
    nexttile(2);
    %yyaxis left
   MI=csvread(['mi_ts/elnino_anom_MI_ver_L' num2str(L) '_lag_' num2str(lag) '_region_' region '.csv']);
    


%
plot(t_sat,MI,'LineWidth',1.5)
%plot(run_ave(MI,12),'LineWidth',1), 
hold on
set(gca,'FontSize',16,'YMinorTick','on')
grid on
grid minor

end

nexttile(1),%legend('spatial lag = 0.25º','spatial lag = 0.50º','spatial lag = 1.0º', ...
    %'spatial lag = 2.0º','sst anom.','location','northwest')
ylabel('$SMI_{WE}$','interpreter','latex')
set(gca,'YLim',[0,2.2],'XTickLabel',{'',''},'TickLabelInterpreter','latex'); grid on; grid minor
% yyaxis right
% set(gca,'YLim',[-5,5])
% fill([t_sat,fliplr(t_sat)],[zeros(size(anom_sat')),-5.*ones(size(anom_sat'))],'k','FaceAlpha',.2,'EdgeColor','none', 'HandleVisibility', 'off')

%yyaxis right
%plot(t_sat(1:end-2),anom_sat(1:end-2),'k-','LineWidth',1)
%annotation('arrow',[.5 .7],[.5 .7],'LineWidth',1.5)
%text(.5, .7,'lag','FontSize',20)
legend('$\delta = 0.25^o$','$\delta = 0.50^o$','$\delta = 1.0^o$', ...
   '$\delta = 2^o$','sst anom.','location','northoutside','Orientation','horizontal','interpreter','latex')

nexttile(2),%yyaxis left
ylabel('$SMI_{NS}$','interpreter','latex')
set(gca,'YLim',[0,2.2],'XTickLabel',{'',''},'TickLabelInterpreter','latex')
grid on; grid minor
% yyaxis right
% set(gca,'YLim',[-5,5])
% fill([t_sat,fliplr(t_sat)],[zeros(size(anom_sat')),-5.*ones(size(anom_sat'))],'k','FaceAlpha',.2,'EdgeColor','none','HandleVisibility','off')

% text(1,2,'a)','FontSize',20)
% text(1,2.2,'b)','FontSize',20)
% text(1,2.3,'c)','FontSize',20)

nexttile(3)
usual = readmatrix(['mi_ts/usual_mi_region_' region '.csv']);
plot(t_sat,usual,'LineWidth',1.5,'Color','k')
ylabel('$SMI_{hist.}$','interpreter','latex'); xlabel('years','interpreter','latex')
%plot(run_ave(MI,12),'LineWidth',1), 
hold on
set(gca,'YLim',[0,2],'FontSize',16,'YMinorTick','on','TickLabelInterpreter','latex')
grid on; grid minor
a_text=text(t_oni(1)-1000,7,'(a)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');
b_text=text(t_oni(1)-1000,4.5,'(b)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');
c_text=text(t_oni(1)-1000,2.1,'(c)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');



a1 = annotation('arrow',[0.5678 0.6021],[0.9091 0.8069],'LineWidth',2);
a2 = annotation('arrow',[0.5654,0.5997],[0.6250,0.5228],'LineWidth',2);


set(gcf,'Position',[306 205 856 528])


%% Comparisons



error=readmatrix('mi_ts/error_region_34.csv');
pears=readmatrix('mi_ts/pearson_region_34.csv');


figure, set(gcf,'Position',[306 205 856 528]);%[600 360 745 422])
tl=tiledlayout(4,1,'TileSpacing','compact');
tl.TileSpacing = 'compact';
%tl.Padding = 'compact';

nexttile, hold on, grid on, box on
plot(t_sat,error,'-','LineWidth',1.5,'Color',[55,126,184]./255)
%plot(t_oni,ONI./4+0.5)
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
ylabel('$AAD(^o)$','Interpreter', 'latex')


nexttile, hold on, grid on, box on
plot(t_sat,pears,'-','LineWidth',1.5,'Color',[55,126,184]./255)
%plot(t_oni,ONI./4+0.5)
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
%legend('ERA5','NOAA OI v2','location','southwest','Interpreter', 'latex')

%xlabel('years','Interpreter', 'latex')
ylabel('$r$','Interpreter', 'latex')

a_text=text(t_oni(1)-1000,2.4,'(a)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');
b_text=text(t_oni(1)-1000,1.1,'(b)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');





error=readmatrix('mi_ts/error_region_gulf.csv');
pears=readmatrix('mi_ts/pearson_region_gulf.csv');

% figure, set(gcf,'Position',[600 360 745 422])
% tl=tiledlayout(2,1,'TileSpacing','compact');
% tl.TileSpacing = 'compact';
%tl.Padding = 'compact';

nexttile, hold on, grid on, box on
plot(t_sat,error,'-','LineWidth',1.5,'Color',[55,126,184]./255)
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
ylabel('$AAD(^o)$','Interpreter', 'latex')

nexttile, hold on, grid on, box on
plot(t_sat,pears,'-','LineWidth',1.5,'Color',[55,126,184]./255)
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1)
%legend('ERA5','NOAA OI v2','location','southwest','Interpreter', 'latex','Orientation','horizontal')

xlabel('years','Interpreter', 'latex')
ylabel('$r$','Interpreter', 'latex')

a_text=text(t_oni(1)-1000,2,'(c)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');
b_text=text(t_oni(1)-1000,1.,'(d)','FontName','Helvetica', 'FontSize',16,'Interpreter', 'latex');


a1 = annotation('arrow',[0.2703 0.2516],[0.9655 0.9016],'LineWidth',2);
a2 = annotation('arrow',[0.4408,0.4221],[0.9692,0.9053],'LineWidth',2);
a3 = annotation('arrow',[0.6803,0.6616],[0.9654,0.9015],'LineWidth',2);
a4 = annotation('arrow',[0.4042,0.3929],[0.5075,0.5758],'LineWidth',2);
a5 = annotation('arrow',[0.4637,0.4525],[0.5018,0.5701],'LineWidth',2);
a6 = annotation('arrow',[0.5467,0.5354],[0.5094,0.5777],'LineWidth',2);





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
    ylabel('hor. SMI')

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
    grid on; box on
   end
% figure(10),legend('ONI','El nino period selected','previous period selected','location','northoutside','orientation','horizontal')
% set(gca,'XTickLabel',years)
%     set(gca,'FontSize',18)
%     set(gcf,'position',[455 420 897 334])
%     ylabel('temperature anomaly')
%     xlabel('years')
% MI=csvread(['../elnino_anom_MI_ver_L4_lag_' num2str(lag) '.csv']);