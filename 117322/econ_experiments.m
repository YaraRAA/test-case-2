clear all;
close all;
data=csvread('econ_experiments_data.csv');
z1 = data(:,1);
n1 = data(:,2);
z2 = data(:,3);
n2 = data(:,4);
theta2 = 2*z2./sqrt(n2);
mu = mean(theta2);
sigma = sqrt(var(theta2));

c = 1.96;
N = quantile(n1,.5)/4;

theta = -1:0.01:1.5;
theta = theta';

p1 = (1/sigma)*normpdf((theta-mu)/sigma).*(normcdf(sqrt(N)*theta-c)+normcdf(-sqrt(N)*theta-c))./(normcdf((sqrt(N)*mu-c)/sqrt(1+N*sigma^2))+normcdf((-sqrt(N)*mu-c)/sqrt(1+N*sigma^2)));
p0 = (1/sigma)*normpdf((theta-mu)/sigma).*(1-normcdf(sqrt(N)*theta-c)-normcdf(-sqrt(N)*theta-c))./(1-normcdf((sqrt(N)*mu-c)/sqrt(1+N*sigma^2))-normcdf((-sqrt(N)*mu-c)/sqrt(1+N*sigma^2)));
prior = (1/sigma)*normpdf((theta-mu)/sigma);
p=plot(theta,prior,'--',theta,p1,'Color',[0 0 0]);
ls=line_fewer_markers(theta,p0,20,'LineStyle','-','Marker','d','MarkerSize',3,'Color','black','LineWidth', .75);
p(1).LineWidth = .75;
p(2).LineWidth = .75;
set(gca,'YTick',[]);
yticks([0]);
yticklabels({' '})
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(gca,'XTick',[-1 -0.5 0 0.5 1 1.5]);
legend('prior','posterior with significance','posterior with no significance','Location','northwest')


%probability of rejection
normcdf((sqrt(N)*mu-c)/sqrt(1+N*sigma^2))+normcdf((-sqrt(N)*mu-c)/sqrt(1+N*sigma^2))


