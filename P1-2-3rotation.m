%% This is the case for Rhea
clear all
dthi=0.02;                   % thita的間隔
dfi=0.02;                    % fi的間隔
ffd=ceil(2*pi/dfi);
z=30                         % 地軸傾角[度]
z=z/180*pi;
dt=60;                       % 每次熱傳導經過時間[sec]
dr=0.01;                     % 每個格點的深度[m]
depth=1;                     % 溫度會變化的深度範圍[m]
Ri=1.436e6;                  % Radius of Iapetus[m]
sigma=5.6704*1e-8;           % Stefan-Boltzman constant [W/s/m^2/K^4]
De=1;                        % Distance from Earth to Sun [AU]
Di=10.11595804;              % Distance from Iapetus to Sun [AU]
Si=1368*((De/Di)^2);         % Solar constant at Iapetus[J/s/m^2]
h=floor(depth/dr)            % 深度軸的格點數
P=200                        % 週期[hour]
P=P*60*60;                   % 週期
n=ceil(P/ffd/dt);            % 每個單位格點的熱傳導次數(每個單位格點經過幾個dt)
% rr=linspace(0,depth,h);

%% Dust mental
lo=1.25e3;       %% Mass density (kg/m^3)
kk=6e-4;         %% thermal conductivity. (J/s/m/K)
Cp=670;           %%specific heat [Joule/kg/K]
Ti=20;            %% Thermal inertia (MKS) (J/(K m^2 s^0.5))
A=0.05
%% Rock
% lo=3;          %% density (g/cm**3)
% kk=6e-3;       %% thermal cond. (W/cm K)
% Cp=670          %%specific heat [Joule/kg.K]
% Ti=1000;        %% Thermal inertia (MKS) (J/K m**2s**0.5)
% A=1.4;
%% constants
dF=Si*(1-A)                 % Flux per unit area
a=zeros(h,h+2);
b=zeros(h,h+2);
R=kk*dt/lo/Cp/((dr)^2)       % R=kk*dt/lo/Cp/((dr)^2)
            for i=1:h
                a(i,i)=-R;
                a(i,i+1)=2+2*R;
                a(i,i+2)=-R;
            end
            a(:,1)=[];
            a(:,h+1)=[];           
            a(1,1)=2+R;
            a(h,h)=2+R;
            for i=1:h
                b(i,i)=R;
                b(i,i+1)=2-2*R;
                b(i,i+2)=R;
            end                
            b(:,1)=[];
            b(:,h+1)=[];      
            b(1,1)=2-R;
            b(h,h)=2-R;
                p=inv(a)*b;             
               
% Solor Flux
td=1;
for thi=-pi/2:dthi:pi/2
    fd=1;
    for fi=0:dfi:2*pi
        F(td,fd)=dF*cos(thi)*sin(fi);
        if(fi>=pi)
            F(td,fd)=0;    %%
        end
        fd=fd+1;
    end
    td=td+1;
end
ftd=td-1;
ffd=fd-1;
clear a b
%% %%%%%%%%%%%%%%%%%%% Main Program%%%%%%%%%%%%%%%%%%%%%%

T=zeros(h,ftd,ffd);  %Temperature(depth,degree of thita,d.. fi)
T(1,:,:)=50;
qT=zeros(h,1);
F=abs(F);                       %計算表面溫度
%% calaulate for many times
hold on;
for time = 1:20
for td=1:ftd
for fd = 1:ffd
    oT=T(2,td,fd);

% 2-divided method for finding roots
T0=0;
T1=2000;
fT0=F(td,fd)+kk/dr*oT-sigma*(T0^4)-kk/dr*T0;
fT1=F(td,fd)+kk/dr*oT-sigma*(T1^4)-kk/dr*T1;
for m=1:30
    if(T1>=10000)
        disp('wrong')
        break;
    end
    if(fT0*fT1>0)
        T1=T1*5;
        fprintf('again in time:%f,fd:%f,T1:%f,fT1:%f,fT0:%f\n',time,fd,T1,fT1,fT0)
        continue;
    end
    T2=(T0+T1)/2;
    fT2=F(td,fd)+kk/dr*oT-sigma*(T2^4)-kk/dr*T2;
    if(fT0*fT2<=0)
        T1=T2;
        fT1=fT2;
    else
        T0=T2;
        fT0=fT2;
    end
    if(abs(fT0-fT1)<=0.01)
        break;
    end
end
for hh=1:n
qT(1)=T1;
qT(:)=p*qT(:);
T(:,td,fd)=qT(:);

end
end
end
if(time==60)
    break;
end
T(:,:,ffd+1)=T(:,:,ffd);
T(:,:,1)=[];
fprintf('%d\n',time);
end
fprintf('done!\n');
plotT=zeros(ftd,ffd);
plotT(:,:)=T(1,:,:);
%%
if (timestart==0)
timestart=time+timestart-2;
end
%%
figure(1)
subplot(1,2,1);
imagesc(plotT(:,:,timestart));
shading interp
hidden off;
colorbar;
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'xtick',gcapoint1);
set(gca,'ytick',gcapoint2);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'-90','0','90'});
title('Temperature Map (P=4.5 days)')
ylabel('Temperature[K]');
xlabel('Local Time');
copy_T=plotT;

%% Clear all except ftd ffd T

clear Cp De Di F P R Ri Si T0 T1 T2 Ti dF depth dfi dr dt dthi fT0 fT1 fT2 fd
clear fi gd grid_time h hh i kk lo m n oT p qT sigma td thi time txt xp z 
clear OP VP VR a b cT
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Constant
 

clear num den_atm num v p_cos p_sin 
K=1.38e-23;      % Boltzmann constant
G=6.67e-11;      % m^3/kg/s^2
Mstar=2.306518e21; % kg
mole=6.02e23;  % moleculer 
distance=10.11595804; % [AU]
dv=1;
den_atm=zeros(ftd,ffd);   % Density of the atmosphere 
lifetime=(distance/0.386)^2*(2.83e6); %[t]
Radii=764.3e3;   % [m] P.C. Thomas et al,,2007
mmass=44e-3;  % Oxygen
ve=sqrt(2*G*Mstar/Radii);          % escape speed [m/s]

%% Probability of sin
qq=1;
n_sin=0;
for i=0:0.01:pi 
    aa=sin(i);
    aa=aa*10000;
    p_sin(n_sin+1:n_sin+1+ceil(aa))=i;
    n_sin=n_sin+ceil(aa);
    qq=qq+1;
end
n_sin=n_sin+1;
%% Probability of cos
qq=1;
n_cos=0;
for i=-pi/2:0.01:pi/2
    aa=cos(i)*1e4;
    p_cos(n_cos+1:n_cos+1+ceil(aa))=i;
    n_cos=n_cos+ceil(aa);
    qq=qq+1;
end
n_cos=n_cos+1;
%%
hold on;
tic;
%figure(3)
flytime=0;
%sphere;
%hold on;
%for thi=-pi/2:0.01:pi/2
%    fi=pi/2;
%   [xx,yy,zz]=sph2cart(-fi+pi/2,thi,1);
%   plot3(xx,yy,zz,'r.');
%   hold on;
%end
for rotation=1:ffd
for nnn=1:10000
thi=p_cos(ceil(rand*n_cos));     % random theta
fi=p_sin(ceil(rand*n_sin));      % random phi
   W=1;
   mm=0;
   alt=0;          % fly time
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   if(fd==0)
       fd=ffd;
   end
   if(td==0)
       td=ftd;
   end
   local_T=plotT(td,fd);
%   fprintf('td=%d,fd=%d,localT=%f',td,fd,local_T);
while(W>=0.01&&local_T>=60)
   if(alt==0)
       den_atm(td,fd)=den_atm(td,fd)-1;
   end
   if(fd>ffd)
       fi=fi-2*pi;
       fd=fd-ffd;
   end
   local_T=plotT(td,fd);
   v=0;
   for mk=1:12
   v=v+rand*12/11.9;
   end
   stdv2=2*K*local_T/mmass*mole;
   v=(v+6)/12*2*sqrt(stdv2);
   if(v-ve>=0)
       fprintf('fly by %f\n',v);
       continue;
   end
   flyd=rand*2*pi;    % 0~pi 飛行方向天頂角
   x=rand*2*pi;
   thi=-thi+pi/2;
   fi=-fi+pi/2;
   d=acos(1-((8*((v/ve)^4)*(sin(x)*cos(x))^2)/(1-4*(v/ve)^2*(1-(v/ve)^2)*sin(x)^2)));
   thi=acos(cos(thi)*cos(d)+sin(thi)*sin(d)*cos(flyd));
   fi=fi+asin(sin(d)*sin(flyd)/sin(thi));
   flytime=v*cos(x)/G/Mstar*(Radii^2);
   if(fi<=pi)
       alt=alt+flytime;
       W=1*exp(-alt/lifetime);
   end
   thi=-thi+pi/2;
   fi=-fi+pi/2;   
   td=ceil((thi+pi/2)/pi*ftd);
   fd=ceil(fi/2/pi*ffd);
   if(fd==0)
       fd=ffd;
   end
   if(td==0)
       fd=ffd;
   end
   if(td<=0)
       if(fd>=ceil(ffd/2))
       fd=fd-ceil(ffd/2);
       fi=fi-pi;
       td=ftd+td;
       thi=thi+pi/2;
       else
       fd=fd+ffd;
       fi=fi+2*pi;
       td=-td;
       thi=-thi;
       end
   end
   if(fd<=0)
       fd=fd+ffd;
       fi=fi+2*pi;
   elseif(fd>ffd)
       fd=fd-ffd;
       fi=fi-2*pi;
   end    
 %  [xx,yy,zz]=sph2cart(-fi-pi/2,thi,1);
 %  plot3(xx,yy,zz,'y.');
 %  hold on;
   mm=mm+1;
end
%fprintf('mm=%d, flytime=%f, localT=%f\n',mm,flytime/3600,local_T);
%plot3(xx,yy,zz,'k*');
den_atm(td,fd)=den_atm(td,fd)+1;
%fprintf('%d\n',nnn);
end
den_atm(:,ffd+1)=den_atm(:,1);
den_atm(:,1)=[];
fprintf('%d\n',rotation);
end
toc;
%%
figure(1)
subplot(1,2,1);imagesc(plotT);colorbar
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'xtick',gcapoint1);
set(gca,'ytick',gcapoint2);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'-90','0','90'});
title('Temperature Map (P=4.5 days)')
ylabel('Latitude[K]');
xlabel('Local Time');
subplot(1,2,2);imagesc(den_atm);
gcapoint1=[1 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
gcapoint2=[1 ceil(ftd/2) ftd];
set(gca,'xtick',gcapoint1);
set(gca,'ytick',gcapoint2);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'-90','0','90'});
title('CO2 distribution [# per 10^4 sputtered Carbon Dioxide]')
ylabel('Latitude[K]');
xlabel('Local Time');
colormap(jet)
colorbar
axis tight
%%
figure(3)
mesh(den_atm);
axis([0 629 0 315 -10 10])
colormap(jet)
colorbar

%%
figure(4)
test=zeros(315,629);
test(:,:)=T(1,:,:);
imagesc(test)
gcapoint1=[0 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
gcapoint2=[0 ceil(ftd/2) ftd];
set(gca,'xtick',[gcapoint1]);
set(gca,'ytick',[gcapoint2]);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'-90','0','-90'});
title('Temperature Map')
xlabel('Local Time');
ylabel('Longitude');
colorbar
hold on;
%%
figure(1)
subplot(1,2,2);
test2=den_atm;
cc=find(ismember(test2,0));
test2(cc)=nan;
surf(test2);
colorbar