%Differential equation dy/dt=-(y+1)(y+3)
clc;
clear ;
%initializing initial value for ode solving
to=0;
tend=2;
yo=-2;
h=0.5;
%initializing of variable for IEM
N=(tend-to)/h;
YE=zeros(N+1,1);
Y=zeros(N+1,1);
t=linspace(to,tend,N+1);
Y(1)=yo;

%implementing implict Euler method 
for i=1:N
    %exact solution of the function
    YE(i)=-3+ 2/(1 + exp(-2*t(i)));
    fi=-((Y(i)+1)*(Y(i)+3));
    Y(i+1)=Y(i)+h*fi;
end
error=Y-YE;
disp(error);
plot(t,Y,t,YE);
%Exact solution 
