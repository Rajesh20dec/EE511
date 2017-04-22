%------------------------------------------------------------------%
%Funtion name:Gamma_Samples
%Input parameters:None
%Output Parameters:None
%Defination:This function Generates 1000 samples of  Gamma 
%distribution using accept reject methods and plot histogram
%and compare it with theoretical PDF
%--------------------------------------------------------------------%
function Gamma_Samples()
clc;
clear all;
No_of_Samples=1000;
Reject_count=0;
for Sample_index=1:No_of_Samples
    while(1)
        u1=rand;
        Y=-5.5*log(u1);
        %below condition is checkin this U2<f(x)/c*g(x)
        %calculation of c and expression for ratio 
        %f(x)/c*g(x) are shown in Project report. 
        if rand<(Y^4.5*exp(-9*Y/11))/((5.5^4.5)*exp(-4.5))
            X(Sample_index)=Y;
            break;
        else
            Reject_count=Reject_count+1;
        end   
    end
end


figure(1)
yyaxis left
hist(X,30);
title('Overlay of Samples and PDF') 
xlabel('x--->')
ylabel('Frequency')
yyaxis right
t=0:0.1:50;
temp=gampdf(t,5.5,1);
plot(t,temp);
ylabel('Gamma  PDF');
legend({'Hist of Random sample','Theoretical PDF of Gamma Distribution'},'FontSize',8)
Acceptance_ratio=No_of_Samples/(No_of_Samples+Reject_count);
sprintf('acceptance rate is =%f',Acceptance_ratio)
end