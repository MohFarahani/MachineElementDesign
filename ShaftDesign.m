function y=ShaftDesign(x)
clc;
clear('all');
p=input('please insert power(Kw): ');
n=input('please insert revolution(rpm): ');
a=input('pressure angel: ');
b=input('other angel: ');
a1=input('please insert Normal Pressure angel for Gear D: ');
ab=input('please insert Length between Bearing A & Gear B:  ');
bc=input('please insert Length between Gear B & Bearing C: ');
cd=input('please insert Length between Bearing C & Gear D: ');
db=input('please insert Diameter of Gear B: ');
dd=input('please insert Diameter of Gear D: ');
t=input('please insert Thickness of the Gear& Bearing: ')
f1=input('please insert Coefficient of diameter between A&B: ');
f2=input('please insert Coefficient of diameter between B&C: ');
f3=input('please insert Coefficient of diameter between C&D: ');
fm=input('please insert Coefficient of diameter for compute the deflection: ');
mb=input('please insert mass of gear b: ');
md=input('please insert mass of gear d: ');
direct=('please insert Direction of Fa: ');
aa=input('please insert a(coefficient of Surface factor(kb=a*Sut^b)): ');
bb=input('please insert b(coefficient of power of Sut in Surface factor(kb=a*Sut^b)): ');
kc=input('please insert Kc (Load modificatian factor): ');
kd=input('please insert Kd(temperature modificatian factor): ');
ke=input('please insert Ke(reliability factor): ');
kff=input('please insert kf (miscellaneous-effects modificatian factor): ');
kfa=input('please insert Kfa(stress concentration for axial stress): ');
kfb=input('please insert Kfb(stress concentration for bending stress): ');
kft=input('please insert Kft(stress concentration for torsional stress): ');
sut=input('please insert Tensile Strength (Mpa): ');
syp=input('please insert Yield Strength(Mpa): ');
fs=input('please insert factor of safety': );
dmax=input('please insert Maximum displacement on Gear(mm):  ');
thetamax=input('please insert Maximum allowable Deflection angel on Bearing(rad): ');
e=input('please insert Modulus of Elasticity(Gpa): ');
w=input('please insert First Critical Speed(rpm): ');
g=9.81; 
e=e*1000;
ka=aa*sut^bb;
sse=0.504*sut;
if sut>1460
    sse=740;
end
a1=a1*pi/180
a=a*pi/180;  %change the radian angle to degree angle
b=b*pi/180;  %change the radian angle to degree angle
T=(p/n)*(30/pi*10^6);  %T=Torsion(N.mm)
ft=2*T/db;
fa=(direct)*ft*tan(b);
fr=ft*tan(a)/cos(b);
f2t=ft*db/dd;          % f2t=F't
f2r=f2t*tan(a1);        % f2r=F'r
ay=-fa;                % ay=bearing force to 'Y' axial
ma=fa*db/2;            % ma=bending to Fa
l=ab+bc+cd;            % l=whole lenght of the shaft
cz=((fr+mb*g)*ab+ma+(md*g-f2r)*l)/(l-cd);   %cz=bearing(C) force to 'z' axial
az=fr+mb*g+md*g-cz-f2r;                     %az=bearing(A) force to 'z' axial
cx=-(ft*ab+f2t*l)/(l-cd);                   %cx=bearing(C) force to 'x' axial
ax=-(ft+f2t+cx);                            %ax=bearing(A) force to 'x' axial
mlb=ab*(((az)^2+(ax)^2)^0.5);               %mlb=Total bending in Right of Gear B
mc=((az*(l-cd)+ma-(mb*g+fr)*bc)^2+(-ax*(l-cd)-ft*(bc))^2)^0.5;    %mc=Total bending in Right of Gear D
% Fatigue method
%
salb=kfb*32*abs(mlb)/(pi*(f1)^3);            
smlb=kfa*4*abs(ay)/(pi*(f1)^2);
smlc=sqrt(3)*kft*16*T/(pi*f2^3);
salc=kfb*32*abs(mc)/(pi*f2^3);
smrc=sqrt(3)*kft*16*T/(pi*f3^3);
sarc=kfb*32*abs(mc)/(pi*f3^3);
q=ka*kc*kd*ke*kff*sse;
%Soderberg Method
%Kf*Sa/Se+Kf*Sm/Syp=1/fs
   % In Left of B 
     f=@(x)(salb/(q*1.24*x^-0.107*x^3)+smlb/(syp*x^3)-1/(fs));
        dsoderbergLB=fsolve(f,10);
          if dsoderbergLB>=51
            f=@(x)salb/(q*1.51*x^-0.157*x^3)+smlb/(syp*x^3)-1/(fs);
            dsoderbergLB=fsolve(f,51);  
          end
   % In Left of C         
     f=@(x)(salc/(q*1.24*x^-0.107*x^3)+smlc/(syp*x^3)-1/(fs));
        dsoderbergLC=fsolve(f,10);
          if dsoderbergLC>=51
            f=@(x)(salc/(q*1.51*x^-0.157*x^3)+smlc/(syp*x^3)-1/(fs));
            dsoderbergLC=fsolve(f,51);  
         end
   % In Right of C    
     f=@(x)(sarc/(q*1.24*x^(-0.107)*x^3)+smrc/(syp*(x^3))-1/(fs));
        dsoderbergRC=fsolve(f,10);
         if dsoderbergRC>=51
            f=@(x)(sarc/(q*1.51*x^-0.157*x^3)+sarc/(syp*x^3)-1/(fs));
            dsoderbergRC=fsolve(f,51); 
         end
   soderberg=[dsoderbergLC dsoderbergRC dsoderbergLB];     
   Dsoderberg=max(soderberg);
   fprintf('Diameter of the shaft(d) according to the Soderberg method is:   %f7.3 mm\n ',Dsoderberg)
%Langer Method
%Kf*Sa/Sy+Kf*Sm/Syp=1/fs
    % In Left of B 
        dlangerLB=((salb+smlb)*fs/syp)^(1/3)/f1;

   % In Left of C         
        dlangerLC=((salc+smlc)*fs/syp)^(1/3)/f2;

   % In Right of C    
        dlangerRC=((sarc+smrc)*fs/syp)^(1/3)/f3;
      
   langer=[dlangerLC dlangerRC dlangerLB];     
   Dlanger=max(langer);
   fprintf('Diameter of the shaft(d) according to the Langer method is:   %f7.3 mm\n ',Dlanger)
   
%Modified Goodman Method
%Kf*Sa/Se+Kf*Sm/Sut=1/fs
   % In Left of B 
     f=@(x)(salb/(q*1.24*x^-0.107*x^3)+smlb/(sut*x^3)-1/(fs));
        dmodifiedgoodmanLB=fsolve(f,10);
          if dmodifiedgoodmanLB>=51
            f=@(x)(salb/(q*1.51*x^-0.157*x^3)+smlb/(sut*x^3)-1/(fs));
           dmodifiedgoodmanLB=fsolve(f,51);  
          end
   % In Left of C         
     f=@(x)(salc/(q*1.24*x^-0.107*x^3)+smlc/(sut*x^3)-1/(fs));
        dmodifiedgoodmanLC=fsolve(f,10);
          if dmodifiedgoodmanLC>=51
            f=@(x)(salc/(q*1.51*x^-0.157*x^3)+smlc/(sut*x^3)-1/(fs));
            dmodifiedgoodmanLC=fsolve(f,51);  
         end
   % In Right of C    
     f=@(x)(sarc/(q*1.24*x^(-0.107)*x^3)+smrc/(sut*(x^3))-1/(fs));
        dmodifiedgoodmanRC=fsolve(f,10);
         if dmodifiedgoodmanRC>=51
            f=@(x)(sarc/(q*1.51*x^-0.157*x^3)+sarc/(sut*x^3)-1/(fs));
            dmodifiedgoodmanRC=fsolve(f,51); 
         end
   modifiedgoodman=[dmodifiedgoodmanLC dmodifiedgoodmanRC dmodifiedgoodmanLB];     
   Ddmodifiedgoodman=max(modifiedgoodman);
   fprintf('Diameter of the shaft(d) according to the Modified Goodman method is:   %f7.3 mm\n ',Ddmodifiedgoodman)
%Gerber Method
%fs*Kf*Sa/Se+(fs*Kf*Sm/Syp)^2=1
   % In Left of B 
     f=@(x)(fs*salb/(q*1.24*x^-0.107*x^3)+(fs*smlb/(sut*x^3))^2-1);
        dgerberLB=fsolve(f,10);
          if dgerberLB>=51
            f=@(x)(fs*salb/(q*1.51*x^-0.157*x^3)+(fs*smlb/(sut*x^3))^2-1);
           dgerberLB=fsolve(f,51);  
          end
   % In Left of C         
     f=@(x)(fs*salc/(q*1.24*x^-0.107*x^3)+(fs*smlc/(sut*x^3))^2-1);
        dgerberLC=fsolve(f,10);
          if dgerberLC>=51
            f=@(x)(fs*salc/(q*1.51*x^-0.157*x^3)+(fs*smlc/(sut*x^3))^2-1);
            dgerberLC=fsolve(f,51);  
         end
   % In Right of C    
     f=@(x)(fs*sarc/(q*1.24*x^(-0.107)*x^3)+(fs*smrc/(sut*(x^3)))^2-1);
        dgerberRC=fsolve(f,10);
         if dgerberRC>=51
            f=@(x)(fs*sarc/(q*1.51*x^-0.157*x^3)+(fs*sarc/(sut*x^3))^2-1);
            dgerberRC=fsolve(f,51); 
         end
         
   gerber=[dgerberLC dgerberRC dgerberLB];     
   Dgerber=max(gerber);
   fprintf('Diameter of the shaft(d) according to the Gerber  method is:   %f7.3 mm\n ',Dgerber)
%ASME-elliptic Method
%(Kf*Sa/Se)^2+(Kf*Sm/Syp)^2=1/fs^2
   % In Left of B 
     f=@(x)((salb/(q*1.24*x^-0.107*x^3))^2+(smlb/(syp*x^3))^2-1/(fs^2));
        dasmeLB=fsolve(f,10);
          if dasmeLB>=51
            f=@(x)((salb/(q*1.51*x^-0.157*x^3))^2+(smlb/(syp*x^3))^2-1/(fs^2));
           dasmeLB=fsolve(f,51);  
          end
   % In Left of C         
     f=@(x)((salc/(q*1.24*x^-0.107*x^3))^2+(smlc/(syp*x^3))^2-1/(fs^2));
        dasmeLC=fsolve(f,10);
          if dasmeLC>=51
            f=@(x)((salc/(q*1.51*x^-0.157*x^3))^2+(smlc/(syp*x^3))^2-1/(fs^2));
            dasmeLC=fsolve(f,51);  
         end
   % In Right of C    
     f=@(x)((sarc/(q*1.24*x^(-0.107)*x^3))^2+(smrc/(syp*(x^3)))^2-1/(fs^2));
        dasmeRC=fsolve(f,10);
         if dasmeRC>=51
            f=@(x)((sarc/(q*1.51*x^-0.157*x^3))^2+(sarc/(syp*x^3))^2-1/(fs^2));
            dasmeRC=fsolve(f,51); 
         end
         
   asme=[dasmeLC dasmeRC dasmeLB];     
   Dasme=max(asme);
   fprintf('Diameter of the shaft(d) according to the ASME-elliptic method is:   %f7.3 mm\n ',Dasme)  
%Rigidity method
   %Deflection in y-z plane
   %EIy''=Az*y-(Fr+Mb*g)<y-ab>^1+Ma*<y-ab>^0+Cz<y-(ab+bc)>^1
   %EIy=Az*y^3/6-(Fr+Mb*g)<y-ab>^3/6+Ma*<y-ab>^2/4+Cz<y-(ab+bc)>^3/6+a1*y+a2
   %a2=0 (Because don't have deflection in bearing A)
   %a1 yield from B.C (Deflection at Bearing C is zero)
   a1=(-az*((l-cd)^3)/6+(fr+mb*g)*(bc^3)/6-ma*(bc^2)/2)/(l-cd);
   dzb=(az*((ab)^3)/6+a1*ab);
   dzd=(az*l^3/6-(fr+mb*g)*(l-ab)^3/6+ma*((l-ab)^2)/2+cz*((cd)^3)/6+a1*l);
   thetaza=a1;
   thetazc=(az*((l-cd)^2)/2-(fr+mb*g)*((bc)^2)/2+ma*(bc)+a1);
   %Deflection in x-y plane
   %EIy''=Ax*y+Ft<y-ab>^1+Cx<y-(ab+bc)>^1
   %EIy=Ax*y^3/6+Ft<y-ab>^3/6+Cx<y-(ab+bc)>^3/6+a2*y+a1
   %a1=0 (Because don't have deflection in bearing A)
   %a2 yield from B.C (Deflection at Bearing C is zero)
   a2=-(ax*((l-cd)^3)/6+ft*((bc)^3)/6)/(l-cd);
   dxb=(ax*((ab)^3)/6+a2*ab);
   dxd=(ax*(l^3)/6+ft*((l-ab)^3)/6+cx*(cd)^3/6+a2*l);
   thetaxa=a2;
   thetaxc=(ax*((l-cd)^2)/2+ft*((bc)^2)/2+a2);
   ddb=(64*sqrt((dzb)^2+(dxb)^2)/(e*dmax*pi))^(1/4)/fm;
   ddd=(64*sqrt((dzd)^2+(dxd)^2)/(e*dmax*pi))^(1/4)/fm;
   thetaA=(64*sqrt((thetaza)^2+(thetaxa)^2)/(e*thetamax*pi))^(1/4)/fm;
   thetaC=(64*sqrt((thetazc)^2+(thetaxc)^2)/(e*thetamax*pi))^(1/4)/fm;
   fprintf('Diameter of the shaft(d) according to the Rigidity (Deflection) method at Gear B:   %f7.3  mm\n',ddb)
   fprintf('Diameter of the shaft(d) according to the Rigidity (Deflection) method at Gear D:   %f7.3  mm\n',ddd)
   fprintf('Diameter of the shaft(d) according to the Rigidity (Deflection Angle)  method at bearing A:   %f7.3  mm\n',thetaA)
   fprintf('Diameter of the shaft(d) according to the Rigidity (Deflection Angle)  method at bearing C:   %f7.3  mm\n',thetaC)
%First Critical Speed
%w^2=(g*(Mb*g*Zb+Md*g*Zd))/(Mb*g*Zb^2+Md*g*Zd^2)
% in this method we only consider weight
  ccz=((mb*ab+md*l)*g)/(l-cd); %Bearing force on C to z axial due to weight
  aaz=(mb+md)*g-ccz;           %Bearing force on A to z axial due to weight
  %EIy''=Azz*y-Mb*g<y-ab>^1+Ccz<y-(ab+bc)>^1
  %EIy=Aaz*y^3/6-Mb*g<y-ab>^3/6+CCz<y-(ab+bc)>^3/6+a3*y+a2
  %a2=0 (Because don't have deflection in bearing A)
  %a3 yield from B.C (Deflection at Bearing C is zero)
  a3=(mb*g*((bc)^3)/6-aaz*(((l-cd)^3)/6))/(l-cd);
  zb=aaz*(ab^3)/6+a3*ab;   %Deflection due to weight on gear A 
  zd=aaz*l^3/6-mb*g*(l-ab)^3/6+ccz*(cd)^3/6+a3*l;  %Deflection due to weight on gear C
  w=2*pi*w/60;
  Dfirstcriticalspeed=(abs(w^2*(mb*g*zb^2+md*g*zd^2)*64/(1000*pi*e*g*(mb*g*zb+md*g*zd))))^(1/4)/fm;
  fprintf('Diameter of the shaft(d) according to the First Critical Speed is:    %f7.3  mm\n',Dfirstcriticalspeed)
%Shear & Bending Force Diagram
% the under function introduce Macoli(Singular) Function 
  ml=@(x,y,z) (((x-y)+abs(x-y))/2).^z.*ceil((atan((((x-y)+abs(x-y))/2))/(pi/2)).^2);
  x=0:0.1:l;
  %Shear Force in Y-Z plane
  %V=Az<y-0>^0-(Fr+Mb*g)<y-ab>^0+Cz<y-(ab+bc)>^0
  vyz=az*ml(x,0,0)-(fr+mb*g)*ml(x,ab,0)+cz*ml(x,ab+bc,0);
  plot(x,vyz)
  title('Shear Force in Y-Z plane')
  xlabel('y(mm)')
  ylabel('V(N)')
  grid on
 %Bending Force in Y-Z plane
 %M=Az<y-0>^1-(Fr+Mb*g)<y-ab>^1+Cz<y-(ab+bc)>^1+Ma<y-ab>^0
  figure
  myz=az*ml(x,0,1)-(fr+mb*g)*ml(x,ab,1)+cz*ml(x,ab+bc,1)+ma*ml(x,ab,0);
  plot(x,myz)
  title('Bending Force in Y-Z plane')
  xlabel('y(mm)')
  ylabel('M(N.mm)')
  grid on
 %Shear Force in X-Y plane
 %V=Ax<y-0>^0+Ft<y-ab>^0+Cx<y-(ab+bc)>^0
  figure
  vxy=ax*ml(x,0,0)+ft*ml(x,ab,0)+cx*ml(x,ab+bc,0);
  plot(x,vxy)
  title('Shear Force in X-Y plane')
  xlabel('y(mm)')
  ylabel('V(N)')
  grid on
  %Bending Force in X-Y plane
  %M=Ax<y-0>^1+Ft<y-ab>^1+Cx<y-(ab+bc)>^1
  figure 
  mxy=ax*ml(x,0,1)+ft*ml(x,ab,1)+cx*ml(x,ab+bc,1);
  plot(x,mxy)
  title('Bending Force in X-Y plane')
  xlabel('y(mm)')
  ylabel('M(N.mm)')
  grid on
%Schematic of the shaft
dlast=[Dsoderberg ddb ddd thetaA thetaC Dfirstcriticalspeed];
d=max(dlast);
yy=[0 ab-t/2 ab-t/2 ab+t/2 ab+t/2 ab+bc ab+bc l-t/2 l-t/2 l+t/2 l+t/2 l-t/2 l-t/2 ab+bc ab+bc ab+t/2 ab+t/2 ab-t/2 ab-t/2 0 0];
z=[f1*d/2 f1*d/2 db/2 db/2 f2*d/2 f2*d/2 f3*d/2 f3*d/2 dd/2 dd/2 -dd/2 -dd/2 -f3*d/2 -f3*d/2 -f2*d/2 -f2*d/2 -db/2 -db/2 -f1*d/2 -f1*d/2 f1*d/2];
figure
plot(yy,z);
axis([-10 l+t -dd dd]);
title('Shematic of the shaft')
xlabel('y(mm)')
ylabel('Z(mm)')
grid on
end 

 