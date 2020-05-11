%Differential equation dy/dt=-(y+1)(y+3)
clc;
clear ;

%initializing initial value for ode solving
to=0;
tend=2;
yo=-2;
h=0.1;

%initializing of variable for MEM
N=(tend-to)/h;
YE=zeros(N+1,1);
Y_EEM=zeros(N+1,1);
Y_MEM=zeros(N+1,1);
Y_IEM=zeros(N+1,1);
t=linspace(to,tend,N+1)';
Y_MEM(1)=yo;
Y_IEM(1)=yo;
Y_EEM(1)=yo;

%implementing implict Euler method 
for i=1:N
    %exact solution of the function
    YE(i)=-3+ 2/(1 + exp(-2*t(i)));
    fi=DIFFE(t(i),Y_EEM(i));
    Y_EEM(i+1)=Y_EEM(i)+h*fi;
end

%implementing Modified Euler method yn+1=yn+h*f(tn+h*0.5,yn+h*0.5*f(yn,tn))
for i=1:N
    %exact solution of the function 
    YE(i)=-3+ 2/(1 + exp(-2*t(i)));
    fi=DIFFE(t(i),Y_MEM(i));
    yn=Y_MEM(i)+h*0.5*(fi);
    Y_MEM(i+1)=Y_MEM(i)+h*(DIFFE(t(i)+h*0.5,yn));
end

%implementing Improved Euler method yn+1=yn+h*0.5*(f(tn, yn)+f(tn+h,yn+h*f(tn,yn))
for i=1:N
    %exact solution of the function 
    fi=DIFFE(t(i),Y_IEM(i));%f(tn,yn)
    yn=Y_IEM(i)+h*fi;%yn+h*f(tn,yn)
    Y_IEM(i+1)=Y_IEM(i)+h*0.5*(fi+DIFFE(t(i)+h,yn));
end
ERR_EEM=abs(Y_EEM-YE);
ERR_MEM=abs(Y_MEM-YE);
ERR_IEM=abs(Y_IEM-YE);

%plotting and finding error
subplot(1,2,1)
plot(t,Y_MEM,t,Y_IEM,t,Y_EEM);
h =legend('eplicit','modified','Improved');
set(h,'Interpreter','none')
title('Implementation of Euler methods')
subplot(1,2,2)
plot(t,ERR_EEM,t,ERR_MEM,t,ERR_IEM);
h =legend('eplicit','modified','Improved');
set(h,'Interpreter','none')
title('Error and Convergence')
%Function to declare differential equation
function x=DIFFE(t,y)
x=-(y+1)*(y+3);
end