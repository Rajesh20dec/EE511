function Proj_P2_Q3()
No_of_samples=input('enter the number of sample required:');
%M_10=1:10;
M_9=0:9;
%Outcomes_M_10=0;
Outcomes_M_9=0;
Outcomes_M_9=datasample(M_9,No_of_samples);
H9= chi2gof(Outcomes_M_9);
hist(Outcomes_M_9);
xlabel('Outcomes');
ylabel('Freqeuncy')
title('Histogram of sampling with replacement');
%Outcomes_M_10=datasample(M_10,No_of_samples);
%H10= chi2gof(Outcomes_M_10);
sprintf('Goodness of fit test at 95 percentage confidence level for samples from 0,1,2...9 is=%d',H9)
Proj_P2_Q3_C(No_of_samples,Outcomes_M_9);
end
function Proj_P2_Q3_C(No_of_samples,Outcomes_M_9)
M_9=0:9;
N0_0f_edges=linspace(1,10,10);
%Outcomes_M_9=datasample(M_9,No_of_samples);
expectedCounts = (No_of_samples* diff(N0_0f_edges));
[h,p,st] = chi2gof(Outcomes_M_9,'edges',N0_0f_edges,'expected',expectedCounts);
sprintf('Goodness of fit test at 95 percentage confidence level for samples from 1,2,3...10 is=%d',h)
end
