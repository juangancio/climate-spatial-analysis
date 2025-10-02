clear all
close all
clc

values = [13.7
5.5
5.5
8
5.5
51.7
19
61
111.5
39.3
95.5
45.8
72.3
305
8
19
61
7.27
13.7
29
14.5
9.9
12.2
66
7.27
10
20
23
41];

h=histogram(values,20);
x = h.BinEdges(1:end-1)+h.BinWidth/2 ; 

hold on
ind = h.Values>0 & x<100 ;
p=polyfit(log(x(ind)),log(h.Values(ind)),1)

x_smooth = 8:300;
plot(x_smooth,x_smooth.^(p(1))*exp(p(2)),'r-','LineWidth',2)

[N,Edges] = histcounts(values,300);
x = Edges(1:end-1)+mean(diff(Edges))/2 ; 
% yyaxis right
% c = cumsum(N);
% plot(x,c./c(end).*100,'k-','LineWidth',1.5)

% c_smooth = cumsum(x_smooth.^(p(1))*exp(p(2)));
% plot(x_smooth,c_smooth./c_smooth(end).*100,'k-')
xlabel('Robustness parameter','interpreter','latex')
ylabel('Number of occurrences','interpreter','latex')
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
grid on