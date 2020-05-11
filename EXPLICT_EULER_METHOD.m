%Differential equation dy/dt=-(y+1)(y+3)
clc;
clear ;

%initializing initial value for ode solving
to=0;
tend=2;
yo=-2;
h=1;

%initializing of variable for EEM
N=(tend-to)/h;
YE=zeros(N+1,1);
Y_EEM=zeros(N+1,1);
t=linspace(to,tend,N+1)';
Y_EEM(1)=yo;

%implementing implict Euler method 
for i=1:N
    %exact solution of the function
    YE(i)=-3+ 2/(1 + exp(-2*t(i)));
    fi=-((Y_EEM(i)+1)*(Y_EEM(i)+3));
    Y_EEM(i+1)=Y_EEM(i)+h*fi;
end
error=abs(YE-Y_EEM);
disp(max(error));
%plotting and finding error
plot(t,Y_EEM,t,YE);

