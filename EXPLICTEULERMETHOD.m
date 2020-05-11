clc;
clear ;
%initializing initial value for ode solving
to=0;
tend=2;
yo=-2;
h=0.1;
N=(tend-to)/h;
y=zeros(N+1,1);
t=linspace(to,tend,N+1);
y(1)=yo;

%implementing implict Euler method 
for i=1:N
    fi=-((y(i)+1)*(y(i)+3));
    y(i+1)=y(i)+h*fi;
end
plot(t,y);
%Exact solution 
