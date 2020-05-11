%Differential equation dy/dt=-(y+1)(y+3)
clc;
clear ;

%initializing initial value for ode solving
to=0;
tend=2;
yo=-2;
h=1;

%initializing of variable for IEM
N=(tend-to)/h;
YE=zeros(N+1,1);
Y_IEM=zeros(N+1,1);
t=linspace(to,tend,N+1)';
Y_IEM(1)=yo;

%implementing Improved Euler method yn+1=yn+h*0.5*(f(tn, yn)+f(tn+h,yn+h*f(tn,yn))
for i=1:N
    %exact solution of the function 
    YE(i)=-3+ 2/(1 + exp(-2*t(i)));
    fi=DIFFE(t(i),Y_IEM(i));%f(tn,yn)
    yn=Y_IEM(i)+h*fi;%yn+h*f(tn,yn)
    Y_IEM(i+1)=Y_IEM(i)+h*0.5*(fi+DIFFE(t(i)+h,yn));
end
%plotting and finding error
error=abs(YE-Y_IEM);
disp(max(error));
plot(t,Y_IEM,t,YE);

%Function to declare differential equation
function x=DIFFE(t,y)
x=-(y+1)*(y+3);
end



