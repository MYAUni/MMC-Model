%Preamble
clc
clear all
clear workspace
%Parameters - Fixed:
mu1=21.05;
gamma=9.12;
p1=0.17;
p2=0.55*10^(-5);
p3=11.9*10^(-7); 
r=0.032;
d0=(1.032*10^5);
mu2=9.12;
k=10^9;
m=1*2395/365;
%Initial conditions
Initial=[ 0 1*10^8 d0/mu2 ]; 
y0 = Initial;
opt = odeset('AbsTol',1e-9,'RelTol',1e-6);
%time of simulation:
tspan1 = linspace(0,365,100); 
%Parameters that vary:
mu1 = 21.05*linspace(0.5,2,11);
a =  100*linspace(0.5,2,11);
%Initial condition 1
for ii = 1:length(mu1)
    for kk = 1:length(a)
%x(1),x(2),x(3) denote the variables M,T,E  respectively.
        odefcn =  @(t,x)[ -mu1(ii)*x(1)+m;
  -x(2)*p1*x(1)/(x(1)+a(kk))+r*x(2)*(1-x(2)/k)-x(2)*(p2*x(3));
  gamma*(p1*x(2)*x(1)/(x(1)+a(kk)))+x(3)*(-mu2)+d0-p3*x(3)*x(2)];
        [t,x(:,:,ii,kk)] = ode45(odefcn, tspan1, y0, opt);
   
    end
end

%Data Extraction:
%This informs us how long our data set is
StoreInside = x(:,:);
%Loop length
Len = length(StoreInside(1,:))/3;
%Storing solutions for T variable:
for ii = 1:Len
    T_1(:,ii) = StoreInside(:,2+3*(ii-1));
    %Skips to next group
end
%Specific times 
for ii = 1:Len
    T1(:,ii) = T_1(length(t),ii);
end
x_contour = a';
y_contour = mu1';
[Y,X] = meshgrid(y_contour,x_contour);
Z1=zeros(length(a),length(mu1));
for i=1:length(a)
  Z1(i,:)=T1(:,(1+(i-1)*length(a)):(i*length(a)));  
end


fig=figure(1);
pcolor(Y,X,Z1); %Another option is to use surf(Y,X,Z1); 
hold on;
shading interp;
cb1 = colorbar;
cb1.LineWidth = 1.5;
ax1=gca;
set(gca,'YScale','log')
xlabel('\mu_{1} [days^{-1}]'), ylabel('a [\muM]');
get(gca,'fontname');  
set(gca,'linewidth', 2,'fontsize',24,'fontname','Helectiva'); % Sets the width of the axis lines, font size, font
set(gca,'TickDir','out');
set(gca,'YScale','log');
h = axes(fig,'visible','off'); 
bottom = min(min(Z1));
top  = max(max(Z1));
clim(ax1,[bottom, top]);                   % Set the same limit to plot
cb1.Position = [0 0 0 0];                   % Reduce the position for cb1
c = colorbar(h,'Position',[0.92 0.23 0.022 0.7]);  % attach colorbar to h
c.LineWidth = 1.5;
c.FontSize = 18;
c.Label.String = 'Cell Count';
clim(cb1.Limits);                          % Set the new colorbar limits




