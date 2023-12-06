
%Preamble
clc
clear all
clear workspace
%Parameters - Fixed:

p1=0.2;
p2=3.7*10^(-6);
d0=1.032*10^5;
mu2=9.12;
a=100;
%%%%%%%%
%Random vectors:
Randomvec1= -4 + (4-(-4)).*rand(5000,1);
Randomvec2= 0.1 + (10-0.1).*rand(5000,1);
%Parameters that vary (scaling a random vector by their model chosen value):
avar = (100.^(Randomvec1))';
mu1var = ((21.05)*Randomvec2)';
rbottom=zeros(1,length(mu1var));
rtop=zeros(1,length(mu1var));
randomr=zeros(1,length(mu1var));


m=zeros(1,length(mu1var));
% Creating the range for parameter r and determining the dose for a
% random chosen value from this range 
for i = 1:length(avar)
    
        a=avar(i);
         mu1=mu1var(i);

          rbottom(i)=p2*d0/mu2;
          rtop(i)=p2*d0/mu2+p1/mu1;
          randomr(i)=rbottom(i)+ (rtop(i)-rbottom(i))*rand(1,1);
          m(i)=mu1*mu1*a*(p2*d0-mu2*randomr(i))/(mu1*mu2*randomr(i)-p2*d0*mu1-mu2*p1);
end

%Scaling for MMC dose

mscaled=365*m;
mls=50;
Molw=334;
conv=10^(6);
mingrams=mscaled/(conv/(Molw*mls));
J=linspace(1,length(Randomvec1),length(Randomvec2));
semilogy(J,mingrams);
xlabel('Index number of the hypothetical patient') 
ylabel('Upper bound value [mg]') 
%hold on



hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = mu1var;
Y =avar;
Z1 = rbottom;
Z2 = rtop;
figure;scatter3(X,Y,Z1,'filled'); hold on;
scatter3(X,Y,Z2,'filled')
line([X; X], [Y; Y], [Z1; Z2], 'Color', 'k');
xlabel('\mu_{1}') 
ylabel('a') 
zlabel('Intervals of parameter r') 





