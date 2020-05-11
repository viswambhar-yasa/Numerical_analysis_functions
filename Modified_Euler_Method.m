%Differential equation dy/dt=-(y+1)(y+3)
clc;
clear ;

%initializing initial value for ode solving
to=0;
tend=2;
yo=-2;
h=1;

%initializing of variable for MEM
N=(tend-to)/h;
YE=zeros(N+1,1);
Y_MEM=zeros(N+1,1);
t=linspace(to,tend,N+1)';
Y_MEM(1)=yo;

%implementing Modified Euler method yn+1=yn+h*f(tn+h*0.5,yn+h*0.5*f(yn,tn))
for i=1:N
    %exact solution of the function 
    YE(i)=-3+ 2/(1 + exp(-2*t(i)));
    fi=DIFFE(t(i),Y_MEM(i));
    yn=Y_MEM(i)+h*0.5*(fi);
    Y_MEM(i+1)=Y_MEM(i)+h*(DIFFE(t(i)+h*0.5,yn));
end
%plotting and finding error
error=abs(YE-Y_MEM);
disp(max(error));
plot(t,Y_MEM,t,YE);

%Function to declare differential equation
function x=DIFFE(t,y)
x=-(y+1)*(y+3);
end



