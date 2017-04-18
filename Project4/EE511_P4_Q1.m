function Proj_4_Q1()
clc
clear all
Proj4_Q1a(10000);
Proj4_Q1b(10000);
Proj4_Q1c(10000);
end

function Proj4_Q1a(No_of_sample)
Summation=0;
for i=1:No_of_sample
    U=rand;
    Summation=(exp((4*U-2)+(4*U-2)^2))*4+Summation;
end
Result=Summation/No_of_sample;
Function=@(x) exp(x+x.*x);
Lower_lim=-2;
Upper_lim=2;
Integral_result=integral(Function,Lower_lim,Upper_lim);
sprintf('Expected value of integration of part a problem =%f and exact value of Integration part a problem =%f',Result,Integral_result)
end
function Proj4_Q1b(No_of_sample)
Summation=0;
for i=1:No_of_sample
    U=rand;
    Summation=Summation-(exp((-1*(1/U-1)*(1/U-1))))*(1/(U.^2));
end
Result=(-2*Summation)/No_of_sample;
Function=@(x) exp(-x.^2);
Lower_lim=0;
Upper_lim=inf;
Integral_result=2*integral(Function,Lower_lim,Upper_lim);
sprintf('Expected value of integration part b problem =%f and exact value of Integration part b problem =%f',Result,Integral_result)
end

function Proj4_Q1c(No_of_sample)
Summation=0;
for i=1:No_of_sample
    U1=rand;
    U2=rand;
    Summation=exp(-1*(U1+U2).*(U1+U2))+Summation;
end
Result=Summation/No_of_sample;
Function=@(x,y) exp(-1*(x+y).*(x+y));
Lower_limX=0;
Upper_limX=1;
Lower_limY=0;
Upper_limY=1;

Integral_result=integral2(Function,Lower_limX,Upper_limX,Lower_limY,Upper_limY);
sprintf('Expected value of integration part c problem=%f and exact value of Integration part c problem=%f',Result,Integral_result)
end


