%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section defines the settings and calls the functions that are to be
% defined later in this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Using matlab solver ode45
close all
clc;
clear;
% Number of steps
N = 5000;
tspan = linspace(0,4*pi,N);
y0 = [1.2; 0; 0; -1.04935750983031990726];
options = odeset('RelTol',1e-3,'AbsTol',1e-3,'Stats','on');

% Settings for figure

figure('position',[100,100,1400,1200],'toolbar','none')
clf

% Calculate solution with ode45 with 'large' tolerance
title('Test 1a: ODE45 with relTol = AbsTol = 0.001.')
[t,y] = ode45(@f,tspan,y0,options);
plotTrajectory(t,y);
%animation(t,y);
pause

title('Test 1a: EULER for 5000 steps')
[t_test,y_test] = eulerEx(@f,tspan,y0); 
plotTrajectory(t_test,y_test);
%animation(t,y);
pause

% Calculate solution with ode45 with smaller tolerance
clf
title('Test 1b: ODE45 with relTol = AbsTol = 0.000001.')
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
[t,y] = ode45(@f,tspan,y0,options); 
plotTrajectory(t,y);
%animation(t,y);
pause

title('Test 1b: EULER for 5000 steps')
[t_test,y_test] = eulerEx(@f,tspan,y0); 
plotTrajectory(t_test,y_test);
%animation(t,y);
pause

% Calculate solution with eulerEx
clf
N = 100000;
tspan = linspace(0,4*pi,N);
title('Test 3: EULER for 10,000 steps')
[tout,yout] = eulerEx(@f,tspan,y0); 
plotTrajectory(tout,yout);
%animation(t,y);
pause

clf
N = 50000;
tspan = linspace(0,4*pi,N);
title('Test 4: EULER for 50,000 steps')
[tout,yout] = eulerEx(@f,tspan,y0); 
plotTrajectory(tout,yout);
%animation(t,y);
pause

clf
N = 5000;
tspan = linspace(0,4*pi,N);
title('Test 5: HEUN for 5,000 steps')
[tout,yout] = heun(@f,tspan,y0); 
plotTrajectory(tout,yout);
%animation(t,y);
pause

clf
N = 20000;
tspan = linspace(0,4*pi,N);
title('Test 6: HEUN for 20,000 steps')
[tout,yout] = heun(@f,tspan,y0); 
plotTrajectory(tout,yout);
animation(t,y);
pause

clf
N = 5000;
tspan = linspace(0,4*pi,N);
y0=[1.3; 0; 0; -1.04935750983031990726];
title('Test 7: ODE45 with relTol = AbsTol = 0.000001. and y0=1.3')
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
[t,y] = ode45(@f,tspan,y0,options); 
plotTrajectory(t,y);
%animation(t,y);
pause

clf
N = 5000;
tspan = linspace(0,4*pi,N);
y0=[1.5; 0; 0; -1.04935750983031990726];
title('Test 8: ODE45 with relTol = AbsTol = 0.000001. and y0=1.5')
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
[t,y] = ode45(@f,tspan,y0,options); 
plotTrajectory(t,y);
%animation(t,y);
pause

clf
N = 5000;
tspan = linspace(0,4*pi,N);
y0=[0.3; 0; 0; -1.04935750983031990726];
title('Test 8: ODE45 with relTol = AbsTol = 0.000001. and y0=0.3')
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
[t,y] = ode45(@f,tspan,y0,options); 
plotTrajectory(t,y);
%animation(t,y);
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END SET UP SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the functions needed 

% function plotTrajectory(~,y) plots earth,moon and the trajectory
function plotTrajectory(~,y)   
    hold on
    plot(0,0,'k.','MarkerSize',40); plot(1,0,'k.','MarkerSize',20);     
    plot(y(:,1),y(:,2),':');   
    axis off   
end

% function animation(t,y) 
function animation(t,y)
    % Animation: simple loop over all time steps
    % Red dot: Satelite, black dots: Earth and Moon.
    % Plot Orbit by blue dots.

    for k = 1:length(t)   
        clf
        hold on
        plot(0,0,'k.','MarkerSize',40); plot(1,0,'k.','MarkerSize',20);     
        plot(y(:,1),y(:,2),':');   
        plot(y(k,1),y(k,2),'r.','MarkerSize',10);
        axis off
        drawnow
    end
end

% Definition of the right-hand side f
function val = f(~,y)    
    % Parameter mu
    mu = 1/82.45; 
     
    %Linear part of f. Matrix.
    A = [0 0  1 0
         0 0  0 1
         1 0  0 2
         0 1 -2 0];
    
    % Non-linear part, splitted into f1 and f2. Column vectors.
    f1 = [zeros(1,size(y,2))
          zeros(1,size(y,2))
          1./((y(1,:)+mu).^2+y(2,:).^2).^1.5.*[y(1,:)+mu
                                      y(2,:) ]];              
       
    f2 = [zeros(1,size(y,2))
          zeros(1,size(y,2))
          1./((y(1,:)-(1-mu)).^2+y(2,:).^2).^1.5.*[y(1,:)-(1-mu)
                                          y(2,:)  ]];
                                      
    % Now f = A*y+(1-mu)*f1 + mu*f2. Result must be column vector.
    val = A*y - (1-mu)*f1 - mu*f2; 
end

%% Define methods as local functions

% explicit Euler
function [tout,yout] = eulerEx(odefun,tspan,y0)
%  eulerEx  Explicite EULER Method with minimalistic interface.
%  [tout,yout] = eulerExplicite(odefun,tspan,y0) odefun function handle
%  f(t,y), tspan vector [t0,te] or vector with given time instances.

    N=length(tspan);
    h=(tspan(N)-tspan(1))/N; 
    yout(:,1)=y0;
    tout(N)=zeros();
    for i=1:N-1
    tout(i)=tspan(i);
    y=yout(:,i);
    t=tout(i);
    k1=odefun(t,y);
    yout(:,i+1)=yout(:,i)+h*(k1);
    end
    yout=yout';
    tout=tout';
end

function [tout,yout] = heun(odefun,tspan,y0,~)
%  heun explicit Heun method with minimalistic interface.
    N=length(tspan);
    h=(tspan(N)-tspan(1))/N; 
    yout(:,1)=y0;
    tout(N)=zeros();
    for i=1:N-1
    tout(i)=tspan(i);
    k1=odefun(tout(i),yout(:,i));
    k2=odefun(tout(i)+(h/3),yout(:,i)+(h/3)*k1);
    k3=odefun(tout(i)+(2*h/3),yout(:,i)+(2*h/3)*k2);
    yout(:,i+1)=yout(:,i)+(h/4)*(k1+3*k3);
    end
    yout=yout';
    tout=tout';
end

