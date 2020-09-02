clear all;
dthi=0.02;                   % thita的間隔
dfi=0.05;                    % fi的間隔
ffd=ceil(2*pi/dfi)
ftd=ceil(pi/dthi);
%z=30                         % 地軸傾角[度]
%z=z/180*pi;
dt=60;                       % 每次熱傳導經過時間[sec]
dr=0.01;                     % 每個格點的深度[m]
depth=1;                     % 溫度會變化的深度範圍[m]
Radii=764.3e3;               % Radius of Rhea[m] P.C. Thomas et al,,2007
sigma=5.6704*1e-8;           % Stefan-Boltzman constant [W/s/m^2/K^4]
De=1;                        % Distance from Earth to Sun [AU]
Di=9.04807635;               % Distance from Saturn to Sun [AU]
Si=1368*((De/Di)^2);         % Solar constant at Iapetus[J/s/m^2]
h=floor(depth/dr)            % 深度軸的格點數
P=4                          % 週期[day]
P=P*24*60*60;                % 週期
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
               
clear a b
%% %%%%%%%%%%%%%%%%%%% Main Program%%%%%%%%%%%%%%%%%%%%%%
T=ones(h,ftd,ffd);  %Temperature(depth,degree of thita,d.. fi)
mole=6.02e23;
T(1,:,:)=50;
qT=zeros(h,1);
timestart=0;
K=1.38066e-23.*mole;
VP=zeros(ftd,ffd); % Vapor Pressure
L=577000; % Latent Heat [J/kg]
VR=zeros(ftd,ffd); % Vaporization rate
mmass=44e-3;             % mass of CO2 [kg/mole]
%% calaulate for many times
hold on;
nturn=1;
ntime=ffd;
for turn=1:nturn
for time = 1:ntime
for td=1:ftd
for fd = 1:ffd
    oT=T(2,td,fd);
    F=abs(Si*(1-A)*cos((td/ftd-0.5)*pi)*sin(fd/ffd*2*pi));
        if(fd>=ceil(ffd/2)),F=0;end
% 2-divided method for finding roots
        T0=0;
        vp=10^(-1367.3/T2+9.9082);
        vr=vp*sqrt(mmass/(2*pi*K*T2));
        fT0=F+kk/dr*oT-sigma*(T0^4)-kk/dr*T0-L*vr;
T1=700;
for m=1:30
%    if(T1>=1000),break;end
%    if(fT0*fT1>0),T1=T1*5;continue;end
    T2=(T0+T1)/2;
    if(T2>138)
        vp=10^(-1275.6/T2+0.00683*T2+8.307);
        vr=vp*sqrt(1/(2*pi*K*T2*mmass));
        fT2=F+kk/dr*oT-sigma*(T2^4)-kk/dr*T2-L*vr;
    else
        vp=10^(-1367.3/T2+9.9082);
        vr=vp*sqrt(1/(2*pi*K*T2*mmass));
        fT2=F+kk/dr*oT-sigma*(T2^4)-kk/dr*T2-L*vr;
    end
    
    if(fT0*fT2<=0),T1=T2;
        else T0=T2;fT0=fT2;end
    if(abs(T0-T1)<=0.01),break;end
end
for hh=1:n
qT(1)=T1;
qT(:)=p*qT(:);
T(:,td,fd)=qT(:);
end
VP(td,fd)=vp;
VR(td,fd)=vr;
end
if(time==ntime&&turn==nturn),break;end
end
T(:,:,ffd+1)=T(:,:,ffd);
T(:,:,1)=[];
fprintf('%d\n',time);
end
end
fprintf('done!\n');
%% Write it on file
plotT=zeros(ftd,ffd);
plotT(:,:)=T(1,:,:);
fid=fopen('RHEA_T_VP_CO2.dat','w');
fprintf(fid,'%d\n',ftd);
fprintf(fid,'%d\n',ffd);
for td=1:ftd
    for fd=1:ffd
        fprintf(fid,'%8.6f\n',plotT(td,fd));
    end
end
fclose(fid);


%%
hold on;
for gd=1:12
    for xp=1:ceil(3*ffd/4)
    plot(xp+ceil(ffd/4),T(1,ceil(ftd/360*30*gd),xp),'b');
    hold on;
    end
    for xp=ceil(3*ffd/4)+1:ffd
    plot(xp-ceil(3*ffd/4),T(1,ceil(ftd/360*30*gd),xp),'b');
    hold on;
    end
    if(gd*20<=ceil(3*ffd/4))
        txt=gd*20+ceil(1/4*ffd);
    else
        txt=gd*20+ceil(3*ffd/4);
    end
    plot(txt,T(1,ceil(ftd/360*30*gd),gd*20),'*b');
    hold on;
    text(txt+5,T(1,ceil(ftd/360*30*gd),gd*20),num2str(gd*30));
end
gcapoint=[0 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
set(gca,'xtick',[gcapoint]);
set(gca,'xticklabel',{'0','6','12','18','24'});
%%
plotT=zeros(ftd,ffd);
plotT(:,:)=T(1,:,:);
surf(plotT(:,:));
shading interp
hidden off;
colorbar;

    

%%
gcapoint1=[0 ceil(ffd/4) ceil(ffd/2) ceil(3*ffd/4) ffd];
gcapoint2=[0 ceil(ftd/2) ftd];
set(gca,'xtick',[gcapoint1]);
set(gca,'ytick',[gcapoint2]);
set(gca,'xticklabel',{'6','12','18','24','6'});
set(gca,'yticklabel',{'-90','0','90'});
title('P=2[hours]')
xlabel('Local Time');
ylabel('Longitude');
zlabel('Temperature[K]');
%%
title('P=200[hours]')
ylabel('Temperature[K]');
xlabel('Local Time');
%% Clear all except ftd ffd T
clear A Cp De Di F P R Ri Si T0 T1 T2 Ti dF depth dfi dr dt dthi fT0 fT1 fT2 fd
clear fi gd grid_time h hh i kk lo m n oT p qT sigma td thi time txt xp z 