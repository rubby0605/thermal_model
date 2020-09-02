disp('The Answer Is,,,')
Radii=1;
 thi=0;
 
    for fi=0:pi/2:3*pi/2
   s(1)=Radii*cos(thi)*sin(fi); % x
   s(2)=Radii*cos(thi)*cos(fi); % y
   s(3)=Radii*sin(thi); 
Ri=sqrt(sum(s(:).^2));
thi=asin(s(3)/Ri);
fi=acot(s(1)/s(2));
   if(s(1)<0),fi=-fi+pi;end
   if(fi<0),fi=fi+pi*3/2;end
   ns(1)=Radii*cos(thi)*sin(fi); % x
   ns(2)=Radii*cos(thi)*cos(fi); % y
   ns(3)=Radii*sin(thi); 
   disp(s);disp(ns);
    end
disp(thi/pi*180)    
%% 
s=[0 1  1]
Ri=sqrt(sum(s(:).^2));
thi=asin(s(3)/Ri);
fi=atan(s(1)/s(2));
   if(s(1)<0),fi=-fi+pi;end
   if(fi<0),fi=fi+pi*3/2;end
      ns(1)=Ri*cos(thi)*sin(fi); % x
      ns(2)=Ri*cos(thi)*cos(fi); % y
      ns(3)=Ri*sin(thi); 
      disp(ns);