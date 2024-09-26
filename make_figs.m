clear all
close all
clc


%% Fig 1

anom=readmatrix('final_ts/elnino_ERA5_monthly_anomaly.csv');
t_oni=datetime(1940,1:length(anom),1);

anom_sat=readmatrix('final_ts/elnino_NOAA_monthly_anomaly.csv');
t_sat=datetime(1981,9:length(anom_sat)+8,1);

anom_std=csvread('final_ts/elnino_ERA5_monthly_std.csv');
anom_sat_std=csvread('final_ts/elnino_NOAA_monthly_std.csv');

anom_gulf=readmatrix('final_ts/gulf_ERA5_monthly_anomaly.csv');
anom_sat_gulf=readmatrix('final_ts/gulf_NOAA_monthly_anomaly.csv');
anom_std_gulf=csvread('final_ts/gulf_ERA5_monthly_std.csv');
anom_sat_std_gulf=csvread('final_ts/gulf_NOAA_monthly_std.csv');


load coastlines.mat
figure

t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding', 'compact');
t2 = tiledlayout(t1,'flow','TileSpacing','Compact','Padding', 'compact');
t3 = tiledlayout(t1,'vertical','TileSpacing','Compact','Padding', 'compact');
t4 = tiledlayout(t1,'vertical','TileSpacing','Compact','Padding', 'compact');


t3.Layout.Tile = 2;
t4.Layout.Tile = 3;
nexttile(t2);
plot(coastlon,coastlat,'k-','LineWidth',1)
set(gca,'YDir','normal')
grid on
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'TickLabelInterpreter','latex')
set(gca,'FontSize',18,'YLim',[-70,70],'XLim',[-180,-20])
hold on



rectangle('Position',[-170,-5,50,10],'FaceColor',[27,158,119]./255) 
rectangle('Position',[-67.5,32.5,22.5,10],'FaceColor',[255,127,0]./255) 

nexttile(t3), hold on

fill([t_oni,fliplr(t_oni)],[anom'-anom_std',fliplr(anom'+anom_std')],'b','FaceAlpha',.4,'EdgeColor','none')
plot(t_oni,anom,'LineWidth',1)

grid on
set(gca,'FontSize',18,'TickLabelInterpreter','latex','XTickLabel',{'','',''})
set(gca,'XLim',[t_sat(1),t_sat(end)])
%figure, hold on
nexttile(t3), hold on
fill([t_sat,fliplr(t_sat)],[anom_sat'-anom_sat_std',fliplr(anom_sat'+anom_sat_std')],'r','FaceAlpha',.4,'EdgeColor','none')
plot(t_sat,anom_sat,'r','LineWidth',1)
grid on,
y=ylabel('SST anomaly  (K)','interpreter','latex','Position',[-1391 4.5000 -1.0000])

set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'XLim',[t_sat(1),t_sat(end)])
set(gcf,'Position',[671 272 934 348])
xlabel('years','Interpreter', 'latex')

a_text=text(t_oni(1)-5000,12.3,'(a)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
b_text=text(t_sat(1)+550,11.5,'(b): El Ni\~{n}o - ERA5','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
c_text=text(t_sat(1)+550,3.25,'(c): El Ni\~{n}o - NOAA OI v2','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');

nexttile(t4), hold on

fill([t_oni,fliplr(t_oni)],[anom_gulf'-anom_std_gulf',fliplr(anom_gulf'+anom_std_gulf')],[255,127,0]./255,'FaceAlpha',.4,'EdgeColor','none')
plot(t_oni,anom_gulf,'LineWidth',1,'Color',[0.8500 0.3250 0.0980])

grid on
set(gca,'FontSize',18,'TickLabelInterpreter','latex','XTickLabel',{'','',''})
set(gca,'XLim',[t_sat(1),t_sat(end)])
%figure, hold on
nexttile(t4), hold on
fill([t_sat,fliplr(t_sat)],[anom_sat_gulf'-anom_sat_std_gulf',fliplr(anom_sat_gulf'+anom_sat_std_gulf')],[255,127,0]./255,'FaceAlpha',.4,'EdgeColor','none')
plot(t_sat,anom_sat_gulf,'LineWidth',1,'Color',[0.8500 0.3250 0.0980])
grid on,%y=ylabel('SST anomaly  (K)','interpreter','latex','Position',[-1.0583e+03 3 -1.0000])

set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'XLim',[t_sat(1),t_sat(end)])
set(gcf,'Position',[671 272 934 348])
xlabel('years','Interpreter', 'latex')
b_text=text(t_sat(1)+550,8.5,'(d): Gulf Stream - ERA5','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
c_text=text(t_sat(1)+200,2.25,'(e): Gulf Stream - NOAA OI v2','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');


%saveas(gcf,'figures/fig1','epsc')

%% Fig 2

sat_ver_lag_1=readmatrix('final_ts/elnino_NOAA_monthly_anom_ver_L4_lag_1.csv');
sat_hor_lag_1=readmatrix('final_ts/elnino_NOAA_monthly_anom_hor_L4_lag_1.csv');

ver_lag_1=readmatrix('final_ts/elnino_ERA5_monthly_anom_ver_L4_lag_1.csv');
hor_lag_1=readmatrix('final_ts/elnino_ERA5_monthly_anom_hor_L4_lag_1.csv');



figure, set(gcf,'Position',[306 205 856 528]);%[600 360 745 422])
tl=tiledlayout(4,1,'TileSpacing','compact');
tl.TileSpacing = 'compact';
%tl.Padding = 'compact';

nexttile, hold on, grid on, box on
plot(t_oni,ver_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
ylabel('$H_{NS}$','Interpreter', 'latex')


nexttile, hold on, grid on, box on
plot(t_oni,hor_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
%legend('ERA5','NOAA OI v2','location','southwest','Interpreter', 'latex')

%xlabel('years','Interpreter', 'latex')
ylabel('$H_{WE}$','Interpreter', 'latex')

a_text=text(t_oni(1)-3500,1.16,'(a)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
b_text=text(t_oni(1)-3500,.8,'(b)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');





sat_ver_lag_1=readmatrix('final_ts/gulf_NOAA_monthly_anom_ver_L4_lag_1.csv');
sat_hor_lag_1=readmatrix('final_ts/gulf_NOAA_monthly_anom_hor_L4_lag_1.csv');

ver_lag_1=readmatrix('final_ts/gulf_ERA5_monthly_anom_ver_L4_lag_1.csv');
hor_lag_1=readmatrix('final_ts/gulf_ERA5_monthly_anom_hor_L4_lag_1.csv');

% figure, set(gcf,'Position',[600 360 745 422])
% tl=tiledlayout(2,1,'TileSpacing','compact');
% tl.TileSpacing = 'compact';
%tl.Padding = 'compact';

nexttile, hold on, grid on, box on
plot(t_oni,ver_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
ylabel('$H_{NS}$','Interpreter', 'latex')

nexttile, hold on, grid on, box on
plot(t_oni,hor_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1)
legend('ERA5','NOAA OI v2','location','southwest','Interpreter', 'latex','Orientation','horizontal')

xlabel('years','Interpreter', 'latex')
ylabel('$H_{WE}$','Interpreter', 'latex')

a_text=text(t_oni(1)-3500,.95,'(c)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
b_text=text(t_oni(1)-3500,.7,'(d)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');

saveas(gcf,'figures/fig2','epsc')

%% Fig 3
sat_ver_lag_8=readmatrix('final_ts/elnino_NOAA_monthly_anom_ver_L4_lag_8.csv');
sat_hor_lag_8=readmatrix('final_ts/elnino_NOAA_monthly_anom_hor_L4_lag_8.csv');

ver_lag_8=readmatrix('final_ts/elnino_ERA5_monthly_anom_ver_L4_lag_8.csv');
hor_lag_8=readmatrix('final_ts/elnino_ERA5_monthly_anom_hor_L4_lag_8.csv');


figure, set(gcf,'Position',[306 205 856 528]);%[600 360 745 422])
tl=tiledlayout(4,1,'TileSpacing','compact');
tl.TileSpacing = 'compact';
%tl.Padding = 'compact';

nexttile, hold on, grid on, box on
plot(t_oni,ver_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
ylabel('$H_{NS}$','Interpreter', 'latex')

varplot=ver_lag_8;
obj=fitlm(1:length(t_oni),varplot);
pval=obj.ModelFitVsNullModel.Pvalue;
lce=obj.Coefficients{2,1};
inter=obj.Coefficients{1,1};
num_t=1:length(t_oni);
plot(t_oni,lce.*num_t+inter,'--','LineWidth',2,'Color',[228,26,28]./255)
str = {['p-value = ' num2str(pval,3) ], ['linear coef. estimate = ' num2str(lce,3)]};
%text(t(40),0.5,str,'FontSize',15,'Color','k','Interpreter','latex')


nexttile, hold on, grid on, box on
plot(t_oni,hor_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18)
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'TickLabelInterpreter','latex','XTickLabel',{'','',''})
%legend('ERA5','NOAA OI v2','location','southwest','Interpreter', 'latex')

% fo = fitoptions('Method','LinearLeastSquares','Robust','LAR');
% varplot=sat_hor_lag_8';
% f=fit(1:length(t_sat),varplot,'poly1',fo)
varplot=hor_lag_8;
obj=fitlm(1:length(t_oni),varplot);
pval=obj.ModelFitVsNullModel.Pvalue;
lce=obj.Coefficients{2,1};
inter=obj.Coefficients{1,1};
num_t=1:length(t_oni);
plot(t_oni,lce.*num_t+inter,'--','LineWidth',2,'Color',[228,26,28]./255)
str = {['p-value = ' num2str(pval,3) ], ['linear coef. estimate = ' num2str(lce,3)]};
%text(t(40),0.83,str,'FontSize',15,'Color','k','Interpreter', 'latex')
%legend('ERA5','NOAA OI v2','location','southwest','Interpreter','latex')

%xlabel('years','Interpreter', 'latex')
ylabel('$H_{WE}$','Interpreter', 'latex')

a_text=text(t_oni(1)-3500,1.35,'(a)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
b_text=text(t_oni(1)-3500,.99,'(b)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');



sat_ver_lag_8=readmatrix('final_ts/gulf_NOAA_monthly_anom_ver_L4_lag_8.csv');
sat_hor_lag_8=readmatrix('final_ts/gulf_NOAA_monthly_anom_hor_L4_lag_8.csv');

ver_lag_8=readmatrix('final_ts/gulf_ERA5_monthly_anom_ver_L4_lag_8.csv');
hor_lag_8=readmatrix('final_ts/gulf_ERA5_monthly_anom_hor_L4_lag_8.csv');


% figure, set(gcf,'Position',[600 360 745 422])
% tl=tiledlayout(2,1,'TileSpacing','compact');
% tl.TileSpacing = 'compact';
%tl.Padding = 'compact';

nexttile, hold on, grid on, box on
plot(t_oni,ver_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1,'XTickLabel',{'','',''})
ylabel('$H_{NS}$','Interpreter', 'latex')

varplot=ver_lag_8;
obj=fitlm(1:length(t_oni),varplot);
pval=obj.ModelFitVsNullModel.Pvalue;
lce=obj.Coefficients{2,1};
inter=obj.Coefficients{1,1};
num_t=1:length(t_oni);
plot(t_oni,lce.*num_t+inter,'--','LineWidth',2,'Color',[228,26,28]./255)
str = {['p-value = ' num2str(pval,3) ], ['linear coef. estimate = ' num2str(lce,3)]};
%text(t(640),0.4,str,'FontSize',15,'Color','k','Interpreter', 'latex')


nexttile, hold on, grid on, box on
plot(t_oni,hor_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'GridColor',[0 0 0],'GridLineWidth',1)
legend('ERA5','NOAA OI v2','location','southwest','orientation','horizontal')

% fo = fitoptions('Method','LinearLeastSquares','Robust','LAR');
% varplot=sat_hor_lag_8';
% f=fit(1:length(t_sat),varplot,'poly1',fo)
varplot=hor_lag_8;
obj=fitlm(1:length(t_oni),varplot);
pval=obj.ModelFitVsNullModel.Pvalue;
lce=obj.Coefficients{2,1};
inter=obj.Coefficients{1,1};
num_t=1:length(t_oni);
plot(t_oni,lce.*num_t+inter,'--','LineWidth',2,'Color',[228,26,28]./255)
str = {['p-value = ' num2str(pval,3) ], ['linear coef. estimate = ' num2str(lce,3)]};
%text(t(640),0.65,str,'FontSize',15,'Color','k','Interpreter', 'latex')
legend('ERA5','NOAA OI v2','location','southwest','Interpreter', 'latex')
xlabel('years','Interpreter', 'latex')
ylabel('$H_{WE}$','Interpreter', 'latex')
a_text=text(t_oni(1)-3500,1.6,'(c)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');
b_text=text(t_oni(1)-3500,.99,'(d)','FontName','Helvetica', 'FontSize',18,'Interpreter', 'latex');

saveas(gcf,'figures/fig3','epsc')

%% Fig. 4

figure
tl=tiledlayout(4,4,'TileSpacing','compact');
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
nexttile(1);
hold on, grid on
sat_ver_lag_1=readmatrix('final_ts/elnino_NOAA_monthly_anom_ver_L4_lag_1.csv');
sat_hor_lag_1=readmatrix('final_ts/elnino_NOAA_monthly_anom_hor_L4_lag_1.csv');

ver_lag_1=readmatrix('final_ts/elnino_ERA5_monthly_anom_ver_L4_lag_1.csv');
hor_lag_1=readmatrix('final_ts/elnino_ERA5_monthly_anom_hor_L4_lag_1.csv');

plot(t_oni,ver_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''}, 'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'YLim',[.3,1],'YTick',[.5,.7,.9])
ylabel('$H_{NS}$','Interpreter', 'latex')
text(t_sat(10),1,'(a)','Interpreter', 'latex','FontSize',18)
title('$lag=0.25^o$','Interpreter', 'latex','FontSize',18)

nexttile(5);
hold on, grid on
plot(t_oni,hor_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''}, 'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
ylabel('$H_{WE}$','Interpreter', 'latex')
text(t_sat(10),1,'(e)','Interpreter', 'latex','FontSize',18)

nexttile(2);
hold on, grid on
sat_ver_lag_2=readmatrix('final_ts/elnino_NOAA_monthly_anom_ver_L4_lag_2.csv');
sat_hor_lag_2=readmatrix('final_ts/elnino_NOAA_monthly_anom_hor_L4_lag_2.csv');

ver_lag_2=readmatrix('final_ts/elnino_ERA5_monthly_anom_ver_L4_lag_2.csv');
hor_lag_2=readmatrix('final_ts/elnino_ERA5_monthly_anom_hor_L4_lag_2.csv');

plot(t_oni,ver_lag_2,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_2,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.3,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(b)','Interpreter', 'latex','FontSize',18)
title('$lag=0.5^o$','Interpreter', 'latex','FontSize',18)

nexttile(6);
hold on, grid on
plot(t_oni,hor_lag_2,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_2,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(f)','Interpreter', 'latex','FontSize',18)

nexttile(3);
hold on, grid on
sat_ver_lag_4=readmatrix('final_ts/elnino_NOAA_monthly_anom_ver_L4_lag_4.csv');
sat_hor_lag_4=readmatrix('final_ts/elnino_NOAA_monthly_anom_hor_L4_lag_4.csv');

ver_lag_4=readmatrix('final_ts/elnino_ERA5_monthly_anom_ver_L4_lag_4.csv');
hor_lag_4=readmatrix('final_ts/elnino_ERA5_monthly_anom_hor_L4_lag_4.csv');

plot(t_oni,ver_lag_4,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_4,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.3,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(c)','Interpreter', 'latex','FontSize',18)
title('$lag=1^o$','Interpreter', 'latex','FontSize',18)

nexttile(7);
hold on, grid on
plot(t_oni,hor_lag_4,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_4,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(g)','Interpreter', 'latex','FontSize',18)

nexttile(4);
hold on, grid on
sat_ver_lag_8=readmatrix('final_ts/elnino_NOAA_monthly_anom_ver_L4_lag_8.csv');
sat_hor_lag_8=readmatrix('final_ts/elnino_NOAA_monthly_anom_hor_L4_lag_8.csv');

ver_lag_8=readmatrix('final_ts/elnino_ERA5_monthly_anom_ver_L4_lag_8.csv');
hor_lag_8=readmatrix('final_ts/elnino_ERA5_monthly_anom_hor_L4_lag_8.csv');

plot(t_oni,ver_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.3,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(d)','Interpreter', 'latex','FontSize',18)
title('$lag=2^o$','Interpreter', 'latex','FontSize',18)

nexttile(8);
hold on, grid on
plot(t_oni,hor_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(h)','Interpreter', 'latex','FontSize',18)

nexttile(9);
hold on, grid on

sat_ver_lag_1=readmatrix('final_ts/gulf_NOAA_monthly_anom_ver_L4_lag_1.csv');
sat_hor_lag_1=readmatrix('final_ts/gulf_NOAA_monthly_anom_hor_L4_lag_1.csv');

ver_lag_1=readmatrix('final_ts/gulf_ERA5_monthly_anom_ver_L4_lag_1.csv');
hor_lag_1=readmatrix('final_ts/gulf_ERA5_monthly_anom_hor_L4_lag_1.csv');

plot(t_oni,ver_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''}, 'FontSize',18,'TickLabelInterpreter','latex')
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
ylabel('$H_{NS}$','Interpreter', 'latex')
text(t_sat(10),1,'(i)','Interpreter', 'latex','FontSize',18)


nexttile(13);
hold on, grid on
plot(t_oni,hor_lag_1,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_1,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'TickLabelInterpreter','latex')
set(gca,'FontSize',18)

set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
%xlabel('years','interpreter','latex')
ylabel('$H_{WE}$','Interpreter', 'latex')
text(t_sat(10),1,'(m)','Interpreter', 'latex','FontSize',18)


nexttile(10);
hold on, grid on
sat_ver_lag_2=readmatrix('final_ts/gulf_NOAA_monthly_anom_ver_L4_lag_2.csv');
sat_hor_lag_2=readmatrix('final_ts/gulf_NOAA_monthly_anom_hor_L4_lag_2.csv');

ver_lag_2=readmatrix('final_ts/gulf_ERA5_monthly_anom_ver_L4_lag_2.csv');
hor_lag_2=readmatrix('final_ts/gulf_ERA5_monthly_anom_hor_L4_lag_2.csv');

plot(t_oni,ver_lag_2,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_2,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(j)','Interpreter', 'latex','FontSize',18)

nexttile(14);
hold on, grid on
plot(t_oni,hor_lag_2,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_2,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'YTickLabel',{'','',''},'TickLabelInterpreter','latex')
set(gca,'FontSize',18)

set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
%xlabel('years','interpreter','latex')
text(t_sat(10),1,'(n)','Interpreter', 'latex','FontSize',18)

nexttile(11);
hold on, grid on
sat_ver_lag_4=readmatrix('final_ts/gulf_NOAA_monthly_anom_ver_L4_lag_4.csv');
sat_hor_lag_4=readmatrix('final_ts/gulf_NOAA_monthly_anom_hor_L4_lag_4.csv');

ver_lag_4=readmatrix('final_ts/gulf_ERA5_monthly_anom_ver_L4_lag_4.csv');
hor_lag_4=readmatrix('final_ts/gulf_ERA5_monthly_anom_hor_L4_lag_4.csv');

plot(t_oni,ver_lag_4,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_4,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(k)','Interpreter', 'latex','FontSize',18)

nexttile(15);
hold on, grid on
plot(t_oni,hor_lag_4,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_4,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'YTickLabel',{'','',''},'TickLabelInterpreter','latex')
set(gca,'FontSize',18)

set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
%xlabel('years','interpreter','latex')
text(t_sat(10),1,'(o)','Interpreter', 'latex','FontSize',18)

nexttile(12);
hold on, grid on
sat_ver_lag_8=readmatrix('final_ts/gulf_NOAA_monthly_anom_ver_L4_lag_8.csv');
sat_hor_lag_8=readmatrix('final_ts/gulf_NOAA_monthly_anom_hor_L4_lag_8.csv');

ver_lag_8=readmatrix('final_ts/gulf_ERA5_monthly_anom_ver_L4_lag_8.csv');
hor_lag_8=readmatrix('final_ts/gulf_ERA5_monthly_anom_hor_L4_lag_8.csv');

plot(t_oni,ver_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_ver_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'XLim',[min(t_sat),max(t_sat)],'XTickLabel',{'','',''},'YTickLabel',{'','',''})
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(l)','Interpreter', 'latex','FontSize',18)

nexttile(16);
hold on, grid on
plot(t_oni,hor_lag_8,'-','LineWidth',1,'Color',[55,126,184]./255)
plot(t_sat,sat_hor_lag_8,'-','LineWidth',1,'Color',[217,95,2]./255)
set(gca,'YLim',[.5,1],'YTick',[.5,.7,.9])
text(t_sat(10),1,'(p)','Interpreter', 'latex','FontSize',18)

set(gcf,'position',[587 232 869 538])



set(gca,'XLim',[min(t_sat),max(t_sat)],'YTickLabel',{'','',''},'TickLabelInterpreter','latex')
set(gca,'FontSize',18)

xlabel(tl,'years','interpreter','latex','FontSize',18)


saveas(gcf,'figures/fig4','epsc')
