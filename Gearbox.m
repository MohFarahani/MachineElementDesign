function y=Gearbox(x)
clc;
clear all;
fs=2;%input('Factor of safety  ');
sut=1770;%input('Tensile Strength  ');
syp=1640;%input('Yield Strength  ');
if sut>1460
    se1=740;
else
    se1=0.504*sut;
end
kd=1;%input('Temperature Factor(Kd) ');
a=4.51;%input('Factor "a" for Surface Factor(Ka) ');
b=-0.265%input('Exponent "b" for Surface Factor(Ka) ');
ka=a*(kd*sut)^b;
ke=0.659;%input('Reliability Factor (Ke) for shaft ');
kf=1;%input('Miscellaneous-Effects Factor (Kf) for  shaft');
kff=1.3;%input('concentration factor (kf)');
phi=20;%input('Pressure angle  ');
phi=phi*pi/180;
phi1=20;%input('Pressure angle for Spur Gear '):
phi1=phi1*pi/180;
sa=15;%input('Helical angle ');
sa=sa*pi/180;
ma=2;%input('module of 5th Gear ');
mb=2;%input('module of 4th Gear ');
mc=1.5;%input('module of 3th Gear ');
md=1.5;%input('module of 2th Gear ');
me=1.5;%input('module of 1th Gear ');
mf=1.5%input('module of Reverse Gear ');
dm=1.5*35;%input('Pinion crownwheel diameter ');
gear5bearing=50; %input('Distance Between Left Top Bearing & 5th gear ');
gear54=50;%input('Distance Between  5th & 4th Gear ');
gear43=50;%input('Distance Between  4th & 3th Gear ');
gear3r=60;%input('Distance Between  3th & Reverse Gear ');
gearR2=60;%input('Distance Between  Reverse Gear & 2th Gear ');
gear21=70; %input('Distance Between 2th & 1th Gear ');
gear1bearing=50;%('Distance Between 1th & Right Top Bearing ');
gear1crownwheel=50;%input('Distance Between 1th Gear & Crownwheel ');
crownwheelbearing=50;%input('Distance Between Right Bottom Bearing & Crownwheel ');
ShaftDistance=114;%input('Distance Between Two shaft');
hp=100;%input('Power(hp)');
n=6000;%input('related Rpm for previous Power ');
%minR1=input('min ratio of Gear H ');
%maxR1=input('max ratio of Gear H ');
%minR2=input('min ratio of Gear F ');
%maxR2=input('max ratio of Gear F ');
%minR3=input('min ratio of Gear E ');
%maxR3=input('max ratio of Gear E ');
%minR4=input('min ratio of Gear D ');
%maxR4=input('max ratio of Gear D ');
%minR5=input('min ratio of Gear B ');
%maxR5=input('max ratio of Gear B ');
mr=mf;%input('Module of Spur Gear for Reverse Gear:  ');
%fprintf('Please insert number of teeth of Reverse Gear On Top Shaft between this two number:  %f7.3 \n',34,40)
Nf=35;%input('');
%fprintf('Please insert number of teeth of Reverse Idlear Gear This between this two number:  %f7.3 \n',27,36)
Nr=30;%input('');
%fprintf('Please insert number of teeth of Reverse Gear On Bottom Shaft between this two number: %f7.3 \n',67,73)
Nl=70;%input('');
thick5=30;%('thick of Gear5');
thick4=35;%('thick of Gear4');
thick3=35;%('thick of Gear3');
thick2=35;%('thick of Gear2');
thick1=35;%('thick of Gear1');
qv5=8;%input('Quality of Kv Factor For Gear5');
qv4=8;%input('Quality of Kv Factor For Gear4');
qv3=8;%input('Quality of Kv Factor For Gear3');
qv2=8;%input('Quality of Kv Factor For Gear2');
qv1=8;%input('Quality of Kv Factor For Gear1');
Y5=0.426;%input('Y cofficient For Lewis For Top Gear 5th ');
Y5b=0.409;%input('Y cofficient For Lewis For Bottom Gear 5th ');
Y4=0.421;%input('Y cofficient For Lewis For Top Gear 4th ');
Y4b=0.415;%input('Y cofficient For Lewis For Bottom Gear 4th ');
Y3=0.43;%input('Y cofficient For Lewis For Top Gear 3th ');
Y3b=0.439;%input('Y cofficient For Lewis For Bottom Gear 3th ');
Y2=0.410;%input('Y cofficient For Lewis For Top Gear 2th ');
Y2b=0.446;%input('Y cofficient For Lewis For Bottom Gear 2th ');
Y1=0.382;%input('Y cofficient For Lewis For Top Gear 1th ');
Y1b=0.5;%input('Y cofficient For Lewis For Bottom Gear 1th ');

life=10^7;%input('Cycle Life For Top Gear');
Rliability=1;%input('Rliability For Gear');
SF=3;%input('Safaty Factor For Bending Strength');
SH=3;%input('Safaty Factor For Contact Strength');
HbGHbP=1;%input=('The Ratio Of Gear Hardness to Pinion Hardness(1 or Between 1.2-1.7)');
EG=191;%input('Modulus Of Elasticity Gear');
EP=2*10^5;%input('Modulus Of Elasticity Pinion');
noGear=0.3;%input('Poissons Ratio For Gear');
noPinion=0.3;%input('Poissons Ratio For Pinion');
JTopGear5=0.62;%input('Bending-Strength Geometry Factor(J) From Figure Top 5th Gear');
JJTopGear5=0.97;%input('Multipliers Of J Factor For Top 5th Gear');
JBottomGear5=0.58;%input('Bending-Strength Geometry Factor(J) From Figure Bottom 5th Gear');
JJBottomGear5=0.98;%input('Multipliers Of J Factor For Bottom 5th Gear');
JTopGear4=0.615;%input('Bending-Strength Geometry Factor(J) From Figure Top 5th Gear');
JJTopGear4=0.985;%input('Multipliers Of J Factor For Top 5th Gear');
JBottomGear4=0.6;%input('Bending-Strength Geometry Factor(J) From Figure Bottom 5th Gear');
JJBottomGear4=0.985;%input('Multipliers Of J Factor For Bottom 5th Gear');
JTopGear3=0.62;%input('Bending-Strength Geometry Factor(J) From Figure Top 5th Gear');
JJTopGear3=1.01;%input('Multipliers Of J Factor For Top 5th Gear');
JBottomGear3=0.64;%input('Bending-Strength Geometry Factor(J) From Figure Bottom 5th Gear');
JJBottomGear3=0.98;%input('Multipliers Of J Factor For Bottom 5th Gear');
JTopGear2=0.58;%input('Bending-Strength Geometry Factor(J) From Figure Top 5th Gear');
JJTopGear2=1.03;%input('Multipliers Of J Factor For Top 5th Gear');
JBottomGear2=0.64;%input('Bending-Strength Geometry Factor(J) From Figure Bottom 5th Gear');
JJBottomGear2=0.97;%input('Multipliers Of J Factor For Bottom 5th Gear');
JTopGear1=0.56;%input('Bending-Strength Geometry Factor(J) From Figure Top 5th Gear');
JJTopGear1=1.01;%input('Multipliers Of J Factor For Top 5th Gear');
JBottomGear1=0.64;%input('Bending-Strength Geometry Factor(J) From Figure Bottom 5th Gear');
JJBottomGear1=0.96;%input('Multipliers Of J Factor For Bottom 5th Gear');
KO=2.25;%input('Overload Factor (KO) For Gear');
cma=[0.14 0.16 0.14 0.14 0.14];

% introduce maculi fanction
ml=@(x,y,z) (((x-y)+abs(x-y))/2).^z.*ceil((atan((((x-y)+abs(x-y))/2))/(pi/2)).^2);
%
%-----------------------------------Design Top Shaft ------------------
%
%Ngmax=(ShaftDistance*2/ma)/(minR5+1);
%Ngmin=ShaftDistance)*2/ma)/(maxR5+1);
%fprintf('Please insert number of teeth of Gear B(not pinion) between this two number %f7.3 \n',Ngmax,Ngmin);
N5b=50 ;%input('Ng=  ');
N5=(ShaftDistance*2/ma)-N5b;
d5=ma*N5;
d5b=ma*N5b;
wt5=(hp*1000*745.7)/(d5*pi*n/60)
wa5=wt5*tan(sa)
wr5=wt5*(tan(phi)/cos(sa))
torqu5=wa5*d5/2

%Nhmax=(ShaftDistance*2/mb)/(minR4+1);
%Nhmin=ShaftDistance)*2/mb)/(maxR4+1);
%fprintf('Please insert number of teeth of Gear B(not pinion) between this two number %f7.3 \n',Nhmax,Nhmin);
N4b=55;%input('Ng=  ');
N4=(ShaftDistance*2/mb)-N4b;
d4=mb*N4;
d4b=mb*N4b;
wt4=(hp*1000*745.7)/(d4*pi*n/60);
wa4=wt4*tan(sa);
wr4=wt4*(tan(phi)/cos(sa));
torqu4=wa4*d4/2;


%Nimax=(ShaftDistance*2/mc)/(minR3+1);
%Nimin=ShaftDistance)*2/mc)/(maxR3+1);
%fprintf('Please insert number of teeth of Gear B(not pinion) between this two number %f7.3 \n',Nimax,Nimin);
N3b=85;%input('Ng=  ');
N3=(ShaftDistance*2/mc)-N3b;
d3=mc*N3;
d3b=mc*N3b;
wt3=(hp*1000*745.7)/(d3*pi*n/60);
wa3=wt3*tan(sa);
wr3=wt3*(tan(phi)/cos(sa));
torqu3=wa3*d3/2;

%Njmax=(ShaftDistance*2/md)/(minR2+1);
%Njmin=ShaftDistance)*2/md)/(maxR2+1);
%fprintf('Please insert number of teeth of Gear B(not pinion) between this two number %f7.3 \n',Njmax,Njmin);
N2b=98;%input('Ng=  ');
N2=(ShaftDistance*2/md)-N2b;
d2=md*N2;
d2b=md*N2b;
wt2=(hp*1000*745.7)/(d2*pi*n/60);
wa2=wt2*tan(sa);
wr2=wt2*(tan(phi)/cos(sa));
torqu2=wa2*d2/2;

%Nkmax=(ShaftDistance*2/me)/(minR1+1);
%Nkmin=ShaftDistance)*2/me)/(maxR1+1);
%fprintf('Please insert number of teeth of Gear B(not pinion) between this two number %f7.3 \n',Nkmax,Nkmin);
N1b=116;%input('Ng=  ');
N1=(ShaftDistance*2/me)-N1b;
d1=me*N1;
d1b=me*N1b;
wt1=(hp*1000*745.7)/(d1*pi*n/60);
wa1=wt1*tan(sa);
wr1=wt1*(tan(phi)/cos(sa));
torqu1=wa1*d1/2;

for j=1:5

se=ka*kd*ke*kf*se1;
wa=[wa5 wa4 wa3 wa2 wa1];
waB=wa;

wt=[wt5 wt4 wt3 wt2 wt1];
wtB=wt;

wr=[wr5 wr4 wr3 wr2 wr1]
wrB=wr;

torqu=[torqu5 torqu4 torqu3 torqu2 torqu1];

Torsion=[wt5*d5/2 wt4*d4/2 wt3*d3/2 wt2*d2/2 wt1*d1/2];

dp=[d5 d4 d3 d2 d1];
dg=[d5b d4b d3b d2b d1b];

LengthOfTopshaft=gear5bearing+gear54+gear43+gear3r+gearR2+gear21+gear1bearing;
x=[gear5bearing gear5bearing+gear54 gear5bearing+gear54+gear43 LengthOfTopshaft-gear1bearing-gear21 LengthOfTopshaft-gear1bearing];

FltbY=(-wa.*dp./2-wr.*(LengthOfTopshaft-x))./LengthOfTopshaft;
FrtbY=-wr-FltbY;
FltbZ=(wt.*(LengthOfTopshaft-x))./LengthOfTopshaft ;   
FrtbZ=wt-FltbZ;
Mleft=FltbY.*x;
Mright=FltbY.*x+torqu;
Mz=FltbZ.*x;
McriticalLeft=((Mleft).^2+(Mz).^2).^0.5
McriticalRight=((Mright).^2+(Mz).^2).^0.5
sBendingAleft=abs(32.*McriticalLeft./pi);
sBendingAright=abs(32.*McriticalRight./pi);
sAxialM=abs(4.*wa./pi);
sTorsionM=abs(16.*Torsion./pi);

sigmaAleft(j)=sBendingAleft(j);
sigmaMleft(j)=sAxialM(j);

% Fatigue method For Left of Gear
  %Soderberg Method
     f=@(x)(sigmaAleft(j)/(se*1.24*x^-0.107*x^3)+sigmaMleft(j)/(syp*x^3)-1/(kff*fs));
        dsoderbergleft(j)=fsolve(f,10);
          if dsoderbergleft(j)>=51
            f=@(x)(sigmaAleft(j)/(se*1.51*x^-0.157*x^3)+sigmaMleft(j)/(syp*x^3)-1/(fs));
            dsoderbergleft(j)=fsolve(f,51);  
          end
 %Langer Method.
        dlangerleft(j)=((kff.*(sigmaAleft(j)+sigmaMleft(j))).*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaAleft(j)./(se*1.24*x^-0.107*x^3)+sigmaMleft(j)./(sut*x^3)-1/(kff*fs));
        dmodifiedgoodmanleft(j)=fsolve(f,10);
          if dmodifiedgoodmanleft(j)>=51
            f=@(x)(sigmaAleft(j)/(se*1.51*x^-0.157*x^3)+sigmaMleft(j)/(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanleft(j)=fsolve(f,51) ;
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaAleft(j)./(se.*1.24.*x^-0.107.*x^3)+(fs.*sigmaMleft(j)./(sut.*x^3)).^2-1);
        dgerberleft(j)=fsolve(f,10);
          if dgerberleft>=51
            f=@(x)(sigmaAleft(j)./(se*1.51.*x^-0.157.*x^3)+(fs.*sigmaMleft(j)./(sut.*x^3))^2-1);
           dgerberleft(j)=fsolve(f,51);
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaAleft(j)./(se.*1.24.*x^-0.107.*x^3)).^2+(sigmaMleft(j)./(syp.*x^3)).^2-1/(fs^2));
        dasmeleft(j)=fsolve(f,10);
        
          if dasmeleft>=51
            f=@(x)((sigmaAleft(j)./(se*1.51*x^-0.157.*x^3)).^2+(sigmaMleft(j)./(syp.*x^3)).^2-1/(fs^2));
           dasmeleft(j)=fsolve(f,51) ;
          end
 %----------- Fatigue method For Right of Gear---------------         
sigmaAright(j)=sBendingAright(j);
sigmaMright(j)=abs(sTorsionM(j));
  %Soderberg Method
     f=@(x)(sigmaAright(j)./(se*1.24*x^-0.107.*x^3)+sigmaMright(j)./(syp*x^3)-1/(kff*fs));
        dsoderbergright(j)=fsolve(f,10);
          if dsoderbergright(j)>=51
            f=@(x)(sigmaAright(j)./(se*1.51.*x^-0.157*x^3)+sigmaMright(j)./(syp.*x^3)-1/(fs));
            dsoderbergright(j)=fsolve(f,51); 
          end
 %Langer Method.
        dlangerright(j)=((kff.*(sigmaAright(j)+sigmaMright(j))).*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaAright(j)./(se*1.24.*x^-0.107.*x^3)+sigmaMright(j)./(sut.*x^3)-1/(kff*fs));
        dmodifiedgoodmanright(j)=fsolve(f,10);
          if dmodifiedgoodmanright(j)>=51
            f=@(x)(sigmaAright(j)./(se*1.51*x^-0.157*x^3)+sigmaMright(j)./(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanright(j)=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaAright(j)./(se*1.24*x^-0.107.*x^3)+(fs.*sigmaMright(j)./(sut.*x^3)).^2-1);
        dgerberright(j)=fsolve(f,10);
          if dgerberright(j)>=51
            f=@(x)(sigmaAright(j)./(se*1.51*x^-0.157.*x^3)+(fs.*sigmaMright(j)./(sut.*x^3)).^2-1);
           dgerberright(j)=fsolve(f,51);  
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaAright(j)./(se*1.24*x^-0.107*x^3)).^2+(sigmaMright(j)./(syp*x^3)).^2-1/(fs^2));
        dasmeright(j)=fsolve(f,10);
        
          if dasmeright(j)>=51
            f=@(x)((sigmaAright(j)./(se*1.51*x^-0.157*x^3)).^2+(sigmaMright(j)./(syp*x^3)).^2-1/(fs^2));
           dasmeright(j)=fsolve(f,51) ; 
          end
%---------------------------------Design Bottom Shaft -------------------------
%
TorsionBottomGear=wt.*dg./2;
wtm=wt.*dg./dm;
wam=wtm.*tan(sa);
wrm=wtm.*(tan(phi)/cos(sa));
TorquG=abs(wa.*dg./2);
TorquM=abs(wam.*dm./2);
LengthOfBottomShaft=gear5bearing+gear54+gear43+gear3r+gearR2+gear21+gear1crownwheel+crownwheelbearing;
p=crownwheelbearing;
FlbbX=abs(wam-wa);
FlbbY=(wr.*(LengthOfBottomShaft-x)-(TorquM+TorquG+wrm.*p))./(LengthOfBottomShaft);
FrbbY=wr-wrm-FlbbY;
FlbbZ=-(wt.*(LengthOfBottomShaft-x)+wtm.*p)/LengthOfBottomShaft;
FrbbZ=-wt-wtm-FlbbZ;
McriticalLeftGear=FlbbY.*x;
McriticalRightGear=FlbbY.*x+TorquG;
McriticalLeftCrownWheel=FlbbY.*(LengthOfBottomShaft-p)+TorquG-wr.*(LengthOfBottomShaft-(x+p));
McriticalRightCrownWheel=McriticalLeftCrownWheel+TorquM;

% Fatigue method For Left of Gear

sBendingAleftGear=abs(32.*McriticalLeftGear./pi);
sBendingArightGear=abs(32.*McriticalRight./pi);
sAxialMleftGear=abs(4.*FlbbX./pi);
sAxialMrightGear=abs(4.*abs(wa+FlbbX)./pi);
sTorsionMrightGear=abs(16.*TorsionBottomGear./pi);

sigmaAleftGear=abs(sBendingAleftGear);
sigmaMleftGear=abs(sAxialMleftGear);
sigmaArightGear=abs(sBendingArightGear);
sigmaMrightGear=((sAxialMrightGear).^2+3.*(sTorsionMrightGear).^2).^0.5;


  %Soderberg Method
     f=@(x)(sigmaAleftGear(j)/(se*1.24*x^-0.107*x^3)+sigmaMleftGear(j)/(syp*x^3)-1/(kff*fs));
        dsoderbergleftGear(j)=fsolve(f,10);
          if dsoderbergleftGear(j)>=51
            f=@(x)(sigmaAleftGear(j)/(se*1.51*x^-0.157*x^3)+sigmaMleftGear(j)/(syp(i)*x^3)-1/(fs));
            dsoderbergleftGear(j)=fsolve(f,51)  ;
          end
 %Langer Method.
        dlangerleftGear(j)=((kff.*(sigmaAleftGear(j)+sigmaMleftGear(j))).*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaAleftGear(j)./(se*1.24*x^-0.107*x^3)+sigmaMleftGear(j)./(sut*x^3)-1/(kff*fs));
        dmodifiedgoodmanleftGear(j)=fsolve(f,10);
          if dmodifiedgoodmanleftGear(j)>=51
            f=@(x)(sigmaAleftGear(j)./(se*1.51*x^-0.157*x^3)+sigmaMleftGear(j)./(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanleftGear(j)=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaAleftGear(j)./(se*1.24.*x^-0.107.*x^3)+(fs.*sigmaMleftGear(j)./(sut*x^3)).^2-1);
        dgerberleftGear(j)=fsolve(f,10);
          if dgerberleftGear(j)>=51
            f=@(x)(sigmaAleftGear(j)./(se*1.51*x^-0.157*x^3)+(fs.*sigmaMleftGear(j)./(sut*x^3)).^2-1);
           dgerberleftGear(j)=fsolve(f,51) ; 
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaAleftGear(j)./(se*1.24*x^-0.107*x^3)).^2+(sigmaMleftGear(j)./(syp*x^3)).^2-1/(fs^2));
        dasmeleftGear(j)=fsolve(f,10);
        
          if dasmeleftGear(j)>=51
            f=@(x)((sigmaAleftGear(j)./(se*1.51*x^-0.157*x^3)).^2+(sigmaMleftGearGear(j)./(syp*x^3)).^2-1/(fs^2));
           dasmeleftGear(j)=fsolve(f,51) ; 
          end
 %----------------- Fatigue method For Right of Gear----------         
  %Soderberg Method
     f=@(x)(sigmaArightGear(j)./(se*1.24*x^-0.107*x^3)+((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(syp)-1/(kff*fs));
        dsoderbergrightGear(j)=fsolve(f,10);
          if dsoderbergrightGear>=51
            f=@(x)(sigmaArightGear(j)./(se*1.51*x^-0.157*x^3)+((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(syp)-1/(kff*fs));
            dsoderbergrightGear(j)=fsolve(f,51);  
          end
 %Langer Method.
        f=@(x)((kff.*(sigmaArightGear(j)./x^3+(kff.*(sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5))-syp/fs);
        dlangerrightGear(j)=fsolve(f,10);
 %Modified Goodman Method
     f=@(x)(sigmaArightGear(j)./(se*1.24*x^-0.107*x^3)+((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(sut)-1/(kff*fs));
        dmodifiedgoodmanrightGear(j)=fsolve(f,10);
          if dmodifiedgoodmanrightGear(j)>=51
            f=@(x)(sigmaArightGear(j)./(se*1.51*x^-0.157*x^3)+((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(sut)-1/(kff*fs));
           dmodifiedgoodmanrightGear(j)=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs*sigmaArightGear(j)./(se*1.24*x^-0.107*x^3)+(fs*((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(sut)).^2-1/kff);
        dgerberrightGear(j)=fsolve(f,10);
          if dgerberrightGear(j)>=51
            f=@(x)(sigmaArightGear(j)./(se*1.51*x^-0.157*x^3)+(fs*((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(sut)).^2-1/kff);
           dgerberrightGear(j)=fsolve(f,51);  
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaArightGear(j)./(se*1.24*x^-0.107*x^3)).^2+(((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(syp)).^2-1/(kff*fs^2));
        dasmerightGear(j)=fsolve(f,10);
        
          if dasmerightGear(j)>=51
            f=@(x)((sigmaArightGear(j)./(se*1.51*x^-0.157*x^3)).^2+(((sAxialMrightGear./x^2).^2+3.*(sTorsionMrightGear./x^3).^2).^0.5./(syp)).^2-1/(kff*fs^2));
           dasmerightGear(j)=fsolve(f,51);  
          end
%------------------Fatigue method For CrownWheel-----------
sBendingAleftCrownWheel=abs(32.*McriticalLeftCrownWheel./pi);
sBendingArightCrownWheel=abs(32.*McriticalRightCrownWheel./pi);
sAxialMleftCrownWheel=sAxialMrightGear;
sTorsionMleftCrownWheel=abs(16.*TorsionBottomGear./pi);

sigmaAleftCrownWheel=abs(sBendingAleftCrownWheel);
sigmaMleftCrownWheel=((sAxialMleftCrownWheel).^2+3.*(sTorsionMleftCrownWheel).^2).^0.5;
sigmaArightCrownWheel=abs(sBendingArightCrownWheel);
sigmaMrightCrownWheel(j)=0;

%----- Fatigue method For Left of CrownWheel---------
  %Soderberg Method
     f=@(x)(sigmaAleftCrownWheel(j)/(se*1.24*x^-0.107*x^3)+((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5/(syp)-1/(kff*fs));
        dsoderbergleftCrownWheel(j)=fsolve(f,10);
          if dsoderbergleftCrownWheel(j)>=51
            f=@(x)(sigmaAleftCrownWheel(j)/(se*1.51*x^-0.157*x^3)+((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5/(syp)-1/(kff*fs));
            dsoderbergleftCrownWheel(j)=fsolve(f,51); 
          end
 %Langer Method.
     f=@(x)(sigmaAleftCrownWheel(j)./(se*1.24*x^-0.107*x^3)+((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5./(syp)-1/(kff*fs));
        dsoderbergleftCrownWheel(j)=fsolve(f,10);
 %Modified Goodman Method
     f=@(x)(sigmaAleftCrownWheel(j)./(se*1.24*x^-0.107*x^3)+((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5./(sut)-1/(kff*fs));
        dmodifiedgoodmanleftCrownWheel(j)=fsolve(f,10);
          if dmodifiedgoodmanleftCrownWheel(j)>=51
            f=@(x)(sigmaAleftCrownWheel(j)/(se*1.51*x^-0.157*x^3)+((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5/(sut)-1/(kff*fs));
           dmodifiedgoodmanleftCrownWheel(j)=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaAleftCrownWheel(j)./(se*1.24*x^-0.107*x^3)+(fs.*((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5./(sut)).^2-1/kff);
        dgerberleftCrownWheel(j)=fsolve(f,10);
          if dgerberleftCrownWheel(j)>=51
            f=@(x)(sigmaAleftCrownWheel(j)./(se*1.51*x^-0.157*x^3)+(fs.*((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5./(sut)).^2-1/kff);
           dgerberleftCrownWheel(j)=fsolve(f,51);  
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaAleftCrownWheel(j)./(se*1.24*x^-0.107*x^3)).^2+(((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5).^2-1/(kff*fs^2));
        dasmeleftCrownWheel(j)=fsolve(f,10);
        
          if dasmeleftCrownWheel(j)>=51;
            f=@(x)((sigmaAleftCrownWheel(j)./(se*1.51*x^-0.157*x^3)).^2+(((sAxialMleftCrownWheel(j)./x^2).^2+3.*(sTorsionMleftCrownWheel(j)./x^3).^2).^0.5./(syp)).^2-1/(kff*fs^2));
           dasmeleftCrownWheel(j)=fsolve(f,51) ; 
          end
 % -----Fatigue method For Right of CrownWheel ----------------        
  %Soderberg Method
     f=@(x)(sigmaArightCrownWheel(j)./(se*1.24*x^-0.107*x^3)+sigmaMrightCrownWheel(j)./(syp*x^3)-1/(kff*fs));
        dsoderbergrightCrownWheel(j)=fsolve(f,10);
          if dsoderbergrightCrownWheel(j)>=51
            f=@(x)(sigmaArightCrownWheel(j)./(se*1.51*x^-0.157*x^3)+sigmaMrightCrownWheel(j)./(syp*x^3)-1/(kff*fs));
            dsoderbergrightCrownWheel(j)=fsolve(f,51) ;
          end
 %Langer Method.
        dlangerrightCrownWheel(j)=((kff.*(sigmaArightCrownWheel(j)+sigmaMrightCrownWheel(j))).*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaArightCrownWheel(j)./(se*1.24*x^-0.107*x^3)+sigmaMrightCrownWheel(j)./(sut*x^3)-1/(kff*fs));
        dmodifiedgoodmanrightCrownWheel(j)=fsolve(f,10);
          if dmodifiedgoodmanrightCrownWheel(j)>=51
            f=@(x)(sigmaArightCrownWheel(j)./(se*1.51*x^-0.157*x^3)+sigmaMrightCrownWheel(j)./(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanrightCrownWheel(j)=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaArightCrownWheel(j)./(se*1.24*x^-0.107*x^3)+(fs.*sigmaMrightCrownWheel(j)./(sut*x^3)).^2-1/kff);
        dgerberrightCrownWheel(j)=fsolve(f,10);
          if dgerberrightCrownWheel(j)>=51
            f=@(x)(fs.*sigmaArightCrownWheel(j)/(se*1.51*x^-0.157*x^3)+(fs.*sigmaMrightCrownWheel(j)/(sut*x^3)).^2-1/kff);
           dgerberrightCrownWheel(j)=fsolve(f,51);  
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaArightCrownWheel(j)./(se*1.24*x^-0.107*x^3)).^2+(sigmaMrightCrownWheel(j)./(syp*x^3)).^2-1/(kff*fs^2));
        dasmerightCrownWheel(j)=fsolve(f,10);
        
          if dasmerightCrownWheel(j)>=51
            f=@(x)((sigmaArightCrownWheel(j)./(se*1.51*x^-0.157*x^3)).^2+(sigmaMrightCrownWheel(j)./(syp*x^3)).^2-1/(kff*fs^2));
           dasmerightCrownWheel(j)=fsolve(f,51);  
          end
end
%          
%---------------------------------Design Shaft Base on Reverse Gear -------------------------      
%----------Top Shaft  
df=mr*Nf;
dr=mr*Nr;
dl=mr*Nl;
rfr=(df+dr)/2;
rlr=(dl+dr)/2;
AlfaFR=acos(((ShaftDistance)^2+rfr^2-rlr^2)/(2*(ShaftDistance)*rfr));
AlfaLR=acos(((ShaftDistance)^2+rlr^2-rfr^2)/(2*(ShaftDistance)*rlr));
wtF=(hp*1000*745.7)/(df*pi*n/60);
wrF=wtF*tan(phi1);
wy=wrF*cos(AlfaFR)+wtF*sin(AlfaFR);
FltbYR=(-wy*(gearR2+gear21+gear1bearing)/LengthOfTopshaft);
FrtbYR=-FltbYR-wy;
wz=wtF*cos(AlfaFR)-wrF*sin(AlfaFR);
FltbZR=wz*(gearR2+gear21+gear1bearing)/LengthOfTopshaft;
FrtbZR=-wz-FltbZR;
McriticalYR=FltbYR*(LengthOfTopshaft-(gearR2+gear21+gear1bearing));
McriticalZR=FltbZR*(LengthOfTopshaft-(gearR2+gear21+gear1bearing));
McriticalR=((McriticalYR).^2+(McriticalZR).^2).^0.5
TorsionR=wt*df/2;
 % Fatigue method         
sigmaAR=32.*McriticalR./pi;
sigmaMR=16.*TorsionR./pi;
  %Soderberg Method
     f=@(x)(sigmaAR/(se*1.24*x^-0.107*x^3)+sigmaMR/(syp*x^3)-1/(kff*fs));
        dsoderbergR(j)=fsolve(f,10)
          if dsoderbergR>=51
            f=@(x)(sigmaAR/(se*1.51*x^-0.157*x^3)+sigmaMR/(syp(i)*x^3)-1/(kff*fs));
            dsoderbergR=fsolve(f,51)
          end
 %Langer Method.
        dlangerR=((kff*(sigmaAR+sigmaMR))*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaAR./(se*1.24*x^-0.107*x^3)+sigmaMR./(sut*x^3)-1/(kff*fs));
        dmodifiedgoodmanR=fsolve(f,10);
          if dmodifiedgoodmanR>=51
            f=@(x)(sigmaAR./(se*1.51*x^-0.157*x^3)+sigmaMR./(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanR=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaAR./(se*1.24*x^-0.107*x^3)+(fs.*sigmaMR./(sut*x^3)).^2-1/kff);
        dgerberR=fsolve(f,10);
          if dgerberR>=51
            f=@(x)(fs.*sigmaAR./(se*1.51*x^-0.157*x^3)+(fs.*sigmaMR./(sut*x^3)).^2-1/kff);
           dgerberR=fsolve(f,51);  
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaAR./(se*1.24*x^-0.107*x^3)).^2+(sigmaMR./(syp*x^3)).^2-1/(kff*fs^2));
        dasmeR=fsolve(f,10);
        
          if dasmeR>=51
            f=@(x)((sigmaAR./(se*1.51*x^-0.157*x^3)).^2+(sigmaMR./(syp*x^3)).^2-1/(kff*fs^2));
           dasmeR=fsolve(f,51);  
          end
% ------------ Bottom Shaft------------------------------      
wtL=wtF;
wrL=wtL*tan(phi1);
wyR=-abs(wrL*cos(AlfaLR)+wtL*sin(AlfaLR));
wtmR=wtL*dl/dm;
wamR=wtmR*tan(sa);
wrmR=wtmR*(tan(phi)/cos(sa));
TorquMR=abs(wamR*dm/2);
FlbbYR=-((wyR*(gearR2+gear21+gear1crownwheel+crownwheelbearing)+wrmR*p+TorquMR)/LengthOfBottomShaft);
FrbbYR=-FlbbYR-wyR-wrmR;
wzR=-wtL*cos(AlfaLR)+wrL*sin(AlfaLR);
FlbbZR=(wtmR*p-wzR*(gearR2+gear21+gear1crownwheel+crownwheelbearing))/LengthOfBottomShaft;
FrbbZR=-wzR-FlbbZR+wtmR;
McriticalReverseYR=FlbbYR*(LengthOfBottomShaft-(gearR2+gear21+gear1crownwheel+crownwheelbearing));
McriticalLeftCrownWheelYR=FlbbYR*(LengthOfBottomShaft-p)+wyR*(gearR2+gear21+gear1crownwheel+crownwheelbearing-p);
TorsionRightReverse=wtL*dl/2


 % Fatigue method   Base on Reverse Gear      
sigmaARightReverse=abs(32*McriticalReverseYR/pi);
sigmaMRightReverse=abs(16*(3*(TorsionRightReverse)^2+(wamR)^2)^0.5/pi);
  %Soderberg Method
     f=@(x)(sigmaARightReverse/(se*1.24*x^-0.107*x^3)+sigmaMRightReverse/(syp*x^3)-1/(kff*fs));
        dsoderbergRightReverse=fsolve(f,10);
          if dsoderbergRightReverse>=51
            f=@(x)(sigmaARightReverse/(se*1.51*x^-0.157*x^3)+sigmaMRightReverse/(syp(i)*x^3)-1/(fs));
            dsoderbergRightReverse=fsolve(f,51); 
          end
 %Langer Method.
        dlangerRightReverse=((kff.*(sigmaARightReverse+sigmaMRightReverse)).*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaARightReverse./(se*1.24*x^-0.107*x^3)+sigmaMRightReverse./(sut*x^3)-1/(kff*fs));
        dmodifiedgoodmanRightReverse=fsolve(f,10);
          if dmodifiedgoodmanRightReverse>=51
            f=@(x)(sigmaARightReverse./(se*1.51*x^-0.157*x^3)+sigmaMRightReverse./(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanRightReverse=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaARightReverse/(se*1.24*x^-0.107*x^3)+(fs.*sigmaMRightReverse/(sut*x^3)).^2-1);
        dgerberRightReverse=fsolve(f,10)
          if dgerberRightReverse>=51
            f=@(x)(fs.*sigmaARightReverse./(se*1.51*x^-0.157*x^3)+(fs.*sigmaMRightReverse/(sut*x^3)).^2-1);
           dgerberRightReverse=fsolve(f,51); 
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaARightReverse./(se*1.24*x^-0.107*x^3)).^2+(sigmaMRightReverse./(syp*x^3)).^2-1/(fs^2));
        dasmeRightReverse=fsolve(f,10);
        
          if dasmeRightReverse>=51
            f=@(x)((sigmaARightReverse./(se*1.51*x^-0.157*x^3)).^2+(sigmaMRightReverse/(syp*x^3)).^2-1/(fs^2));
           dasmeRightReverse=fsolve(f,51);  
          end

 % Fatigue method    Base on CrownWheel     
sigmaALeftCrownWheel=abs(32.*McriticalLeftCrownWheel/pi);
sigmaMLeftCrownWheel=abs(16.*(3.*(TorsionRightReverse).^2+(wamR).^2).^0.5/pi);
  %Soderberg Method
     f=@(x)(sigmaALeftCrownWheel./(se*1.24*x^-0.107*x^3)+sigmaMLeftCrownWheel./(syp*x.^3)-1/(kff*fs));
        dsoderbergLeftCrownWheel=fsolve(f,10);
          if dsoderbergLeftCrownWheel>=51
            f=@(x)(sigmaALeftCrownWheel./(se*1.51*x^-0.157*x^3)+sigmaMLeftCrownWheel./(syp*x^3)-1/(fs));
            dsoderbergLeftCrownWheel=fsolve(f,51); 
          end
 %Langer Method.
        dlangerLeftCrownWheel=((kff*(sigmaALeftCrownWheel+sigmaMLeftCrownWheel)).*fs/syp).^(1/3);
 %Modified Goodman Method
     f=@(x)(sigmaALeftCrownWheel./(se*1.24*x^-0.107*x^3)+sigmaMLeftCrownWheel./(sut*x^3)-1/(kff*fs));
        dmodifiedgoodmanLeftCrownWheel=fsolve(f,10);
          if dmodifiedgoodmanLeftCrownWheel>=51
            f=@(x)(sigmaALeftCrownWheel./(se*1.51*x^-0.157*x^3)+sigmaMLeftCrownWheel./(sut*x^3)-1/(kff*fs));
           dmodifiedgoodmanLeftCrownWheel=fsolve(f,51);  
          end       

 %Gerber Method
     f=@(x)(fs.*sigmaALeftCrownWheel./(se*1.24*x^-0.107*x^3)+(fs.*sigmaMLeftCrownWheel./(sut*x^3))^2-1);
        dgerberLeftCrownWheel=fsolve(f,10);
          if dgerberLeftCrownWheel>=51
            f=@(x)(sigmaALeftCrownWheel./(se*1.51*x^-0.157*x^3)+(fs.*sigmaMLeftCrownWheel./(sut*x^3))^2-1);
           dgerberLeftCrownWheel=fsolve(f,51);  
          end
 %ASME-elliptic Method
   
     f=@(x)((sigmaALeftCrownWheel./(se*1.24*x^-0.107*x^3)).^2+(sigmaMLeftCrownWheel./(syp*x^3)).^2-1/(fs^2));
        dasmeLeftCrownWheel=fsolve(f,10);
        
          if dasmeLeftCrownWheel>=51
            f=@(x)((sigmaALeftCrownWheel./(se*1.51*x^-0.157*x^3)).^2+(sigmaMLeftCrownWheel./(syp*x^3)).^2-1/(fs^2));
           dasmeLeftCrownWheel=fsolve(f,51);  
          end




          
 
WT=[wt5 wt4 wt3 wt2 wt1];
thick=[thick5 thick4 thick3 thick2 thick1];
dTopGear=[d5 d4 d3 d2 d1];
dBottomGear=[d5b d4b d3b d2b d1b];
TeethTopGear=[N5 N4 N3 N2 N1];
TeethBottomGear=[N5b N4b N3b N2b N1b];
QV=[qv5 qv4 qv3 qv2 qv1];
B=0.25.*(12-QV).^2/3;
A=50+56.*(1-B);
VTopGear=n*pi.*dTopGear*0.001;
KVTopGear=((A+(200.*VTopGear).^0.5)./A).^B;
VBottomGear=n*pi.*dBottomGear*0.001;
KVBottomGear=((A+(200.*VBottomGear).^0.5)./A).^B;
ModuleNormal=[ma mb mc md me].*cos(sa);
YLewisTopGear=[Y5 Y4 Y3 Y2 Y1];
YLewisBottomGear=[Y5b Y4b Y3b Y2b Y1b];
KSTopGear=0.8433.*(thick.*ModuleNormal.*(YLewisTopGear).^0.5);
KSBottomGear=0.8433.*(thick.*ModuleNormal.*(YLewisBottomGear).^0.5);
if thick<=25
CpfTopGear=thick/(10.*dTopGear)-0.025;
end

if 25<thick<=432
CpfTopGear=thick/(10.*dTopGear)-0.0375+0.5*0.001.*thick;
end
if 432<thick<=1000
CpfTopGear=thick/(10.*dTopGear)-0.1109+0.815*0.001.*thick-0.3534*10^(-6).*(thick).^2 ;
end
if thick<=25
CpfBottomGear=thick/(10.*dBottomGear)-0.025;
end

if 25<thick<=432
CpfBottomGear=thick/(10.*dBottomGear)-0.0375+0.5*0.001.*thick;
end
if 432<thick<=1000
CpfBottomGear=thick/(10.*dBottomGear)-0.1109+0.815*0.001.*thick-0.3534*10^(-6).*(thick).^2 ;
end
KMTopGear=1+(CpfTopGear+cma);
KMBottomGear=1+(CpfBottomGear+cma);
Pn=pi.*[ma mb mc md me].*cos(sa);
PN=Pn.*cos(phi);
aa=[ma mb mc md me].*cos(sa);
rbpa=(dTopGear+aa).^2-dTopGear.*cos(tan(phi)/cos(sa))/2;
rbga=(dBottomGear+aa).^2-dBottomGear.*cos(tan(phi)/cos(sa));
rpg=0.5.*(dBottomGear+dTopGear).*sin(tan(phi)/cos(sa));
if rbpa<rpg
    rbpa=rbg
end
if rbga<rpg
    rbga=rbg
end
Z=(rbpa).^0.5+(rbga).^0.5-rpg;
MN=PN/(0.95.*Z);
for qq=1:5
MG(qq)=TeethBottomGear(qq)./TeethTopGear(qq);
if MG(qq)<1
    MG=MG(qq).^(-1);
end
end
I=sin(2*tan(phi)/cos(sa))./MN.*MG./(MG+1);
YNTopGear=1.3558*(life)^-0.0178;
YNBottomGear=1.3558.*(dTopGear.*life./dBottomGear).^-0.0178;
ZNTopGear=1.4488*(life)^-0.023;
ZNBottomGear=1.4488.*(dTopGear.*life./dBottomGear).^-0.023;
if Rliability<0.99
   KR=0.658-0.0759*log(1-Rliability);
end
if 0.99<=Rliability<=0.99999
KR=0.5-0.109*log(1-Rliability);

end
if Rliability==1
    KR=1;
end
if HbGHbP==1
    CH=1;
else 
    CHGear=1+(8.98*0.001*HbGHbP-8.29*0.001).*(MG-1);
end
CP=(1/(pi*((1-(noPinion)^2)/EP+(1-(noGear)^2)/EG)))^0.5;
JFactorTopGear=[JTopGear5 JTopGear4 JTopGear3 JTopGear2 JTopGear1].*[JJTopGear5 JJTopGear4 JJTopGear3 JJTopGear2 JJTopGear1];
JFactorBottomGear=[JBottomGear5 JBottomGear4 JBottomGear3 JBottomGear2 JBottomGear1].*[JJBottomGear5 JJBottomGear4 JJBottomGear3 JJBottomGear2 JJBottomGear1];
mtt=[ma mb mc md me];
SigmaBendingTopGear=WT.*KO.*KVTopGear.*KSTopGear.*KMTopGear./(thick.*JFactorTopGear.*mtt);
SigmaBendingBottomGear=WT.*KO.*KVBottomGear.*KSBottomGear.*KMBottomGear./(thick.*JFactorBottomGear.*mtt);
for qq=1:5
    if dBottomGear(qq)<dTopGear(qq);
        dPinion(qq)=dBottomGear(qq);
    else
       dPinion(qq)=dTopGear(qq);
    end
end

SigmaContactTopGear=CP.*(WT.*KO.*KVTopGear.*KSTopGear.*KMTopGear./(dPinion.*thick.*I)).^2
SigmaContactBottomGear=CP.*(WT.*KO.*KVBottomGear.*KSBottomGear.*KMBottomGear./(dPinion.*thick.*I)).^2

SigmaFPpinion=SigmaBendingTopGear.*SF.*KR./YNTopGear
SigmaFPgear=SigmaBendingBottomGear.*SF.*KR./YNBottomGear

SigmaHPpinion=SigmaContactTopGear.*SF.*KR./(CH.*ZNTopGear)
SigmaHPgear=SigmaContactBottomGear.*SF.*KR./(CH.*ZNBottomGear)
         
          
          







