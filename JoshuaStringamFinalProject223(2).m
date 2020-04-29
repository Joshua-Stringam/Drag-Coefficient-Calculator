
clear; clc
%Prompt Variables for 1st guess

Velocity= 3100 * 0.3048; %m/s
Cd= .5; %1st guess 
m= 55*0.000064799;%Projectile Mass in kg (Grains conversion)
d=0.223; %Calibre in fractions of an inch
A= (pi*(d/2)^2)*0.00064516;%Cross sectional area of projectile
g= 9.807; %m/s^2
PT= 102000; %Total pressure in Pa or (N/m^2)
T= 273.15+25; %Temperature in Kelvin
H= 0.49; %Relative Humidity as fraction
%Imperacle data for saturation pressure of water given temperature
PS=0.16004123*(T^3)-139.71796235*(T^2)+40769.2199177*T-3973384.2924139;
PW=PS*H;
PA=PT-PW;
%Effective Density of air at temp, pressure and humidity
ro=(PA/(287.05*T))+(PW/(461.495*T));
%Desired number of iterations to fit CD and angle
PMAX=100;
P=PMAX;%Specify desired number of Iterations
            %Prompt values for data match
            xm=[0, 100, 200, 300, 400, 500];
            ym=[-1.5, 1.7, 0, -8.2, -25.3, -55.1];
            %Estimate angle based on first two data points
            Angle= atan((ym(2)-ym(1))*0.02777778/(xm(2)-xm(1)));
%Allocate space for variables to estimate derivatives for newton's Method
SESQ3=0;
SESQ2=0;
SESQ1=0;

for IT=1:5;
for p=1:5

%Drag Equation
%Fd=0.5*ro*u^2*Cd*A

%%for p=1:P
    %Time Loop
    tfinal=1; %Calculation Duration Seconds
    h=0.0002; %time interval seconds
    N=tfinal/h;

    %Preallocate space for x,y,v dx, dy, dv
    x= zeros(1,N);
    y= zeros(1,N);
    v= zeros(1,N);
    a= zeros(1,N);
    vx= zeros(1,N);
    vy= zeros(1,N);
    dx= zeros(1,N);
    dy= zeros(1,N);
    dv= zeros(1,N);

    %Initial values
    k= 0.5*ro*Cd*A*(1/m);%Deceleration constant, includes Cd, A ro, mass

            %Loop through Data to aquire velocity values over time
            a(1)= Angle;
            v(1)= Velocity(1);
            vx(1)=cos(a(1))*v(1);
            vy(1)=sin(a(1))*v(1);
            dv(1)=(v(1)^2)*k;
            dx(1)=-cos(a(1))*dv(1);
            dy(1)=(-sin(a(1))*dv(1))-9.807;
            a(1)=atan(vy(1)/vx(1));

            for n=1:N-1
            vx(n+1)=vx(n)+(h*dx(n));
            vy(n+1)=vy(n)+(h*dy(n));
            v(n+1)=((vx(n+1)^2)+(vy(n+1)^2))^0.5;
            dv(n+1)=(v(n+1)^2)*k;
            a(n+1)=atan(vy(n+1)/vx(n+1));
            dx(n+1)=-cos(a(n+1))*dv(n+1);
            dy(n+1)=-sin(a(n+1))*dv(n+1)-9.807;
            end

            %Integration of Velocity data (Simpsons Rule) to obtain coordinate data
            x(1)=0;
            y(1)=-0.0381;% Barrel below scope
            for n=1:N-1
            x(n+1)=(((2/3)*(x(n)+vx(n)*h))+((1/3)*(x(n)+(h/2)*(vx(n)+vx(n+1)))));
            y(n+1)=(((2/3)*(y(n)+vy(n)*h))+((1/3)*(y(n)+(h/2)*(vy(n)+vy(n+1)))));
            end
            %%Convert values from meters to inches and yards
            y=y.*39.37;
            x=x.*1.0936;
            if IT==1
            if p==4;
            plot(xm,ym,'o')
            hold on
            plot (x,y,'--');
            hold on
            end
            end
            if IT==2
            if p==4;
            plot (x,y,'--');
            hold on
            end 
            end
            if IT==3
            if p==4;
            plot (x,y);
            hold on
            end
            end
            if IT==5
            if p==4;
            plot (x,y);
            hold on
            end
            end
            L=length(xm);
            error=zeros(1,L);
    for l=1:L
        result(l)=interp1(x,y,xm(l));
        error(l)=ym(l)-result(l);
        ESQ(l)=error(l)^2;
    end
  %Newton's method for calculating drag coefficcient
  DCD=0.002;
  if p==1
  SESQ2=sum(ESQ);
  Cd=Cd-DCD;
  end
  if p==2
  SESQ1=sum(ESQ);
  Cd=Cd+(2*DCD);
  end
  if p==3
  SESQ3=sum(ESQ);
  FP=(SESQ3-SESQ1)/(2*DCD);
  FPP=(SESQ3-2*SESQ2+SESQ1)/(DCD^2);
  Cmod=FP/FPP;
  Cd=Cd-DCD-Cmod;
  end
  
   if p==4;
   %Readjust angle to line up with second data point
   Angle=Angle+((ym(2)-(result(2)))*0.02777778/(xm(2)-xm(1)));
   end
   end
end
for l=1:L
        result(l)=interp1(x,y,xm(l));
        Velocity(l)=interp1(x,v,xm(l));
        Energy(l)= (Velocity(l)^2)*m*(0.305^2)*3.962862;    
end
%%end
xlabel('Yards Distance')
ylabel('Inches Vertical Displacement')
title('Drop Chart Data 0.223 calibre, 55 Grain 3100ft/s')
legend('Experimental Data','1st Guess','2nd guess','3rd guess','5th guess')

