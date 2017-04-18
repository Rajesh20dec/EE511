function Proj_P4_Q2()
clc
clear all
x=[];
sprintf('Below printinf for 10 samples')
[x Theoretical_Prob_dist]=Emperical_distribution(10);
figure(1)
cdfplot(x)
title('Empirical Vs Theoretical Prob.Distribution');
xlabel('x')
ylabel('F(x)')
hold on
grid on
plot(1:10,Theoretical_Prob_dist,'r','linewidth',2);
legend('Empirical dist','Theoretical Prob.Distribution')
hold off


sprintf('Below printinf for 100 samples')
[x Theoretical_Prob_dist]=Emperical_distribution(100);
figure(2)
cdfplot(x)
title('Empirical Vs Theoretical Prob.Distribution');
xlabel('x')
ylabel('F(x)')
hold on
grid on
plot(1:100,Theoretical_Prob_dist,'r','linewidth',2);
legend('Empirical dist','Theoretical Prob.Distribution')
hold off


sprintf('Below printinf for 1000 samples')
[x Theoretical_Prob_dist]=Emperical_distribution(1000);
figure(3)
title('Empirical Vs Theoretical Prob.Distribution');
xlabel('x')
ylabel('F(x)')
cdfplot(x)
hold on
grid on
plot(1:1000,Theoretical_Prob_dist,'r','linewidth',2);
legend('Empirical dist','Theoretical Prob.Distribution')
hold off

end


function[x Theoretical_Prob_dist]= Emperical_distribution(No_of_samples)
X=[];
%No_of_samples=10;
for i=1:No_of_samples
Z=randn(1,4);
X(i)=sum(Z.*Z);
end
x = sort(X);
%below line calculate the emperical disribution from Sample of random
%variable
[Emp_Prob_Dist Samples]=ecdf(x);
Emp_Prob_Dist=Emp_Prob_Dist.';
Samples=Samples.';
Emp_Prob_Dist=Emp_Prob_Dist(2:No_of_samples+1);
Samples=Samples(2:No_of_samples+1);
Theoretical_Prob_dist=chi2cdf(x,4);
Lower_bound=max(Emp_Prob_Dist-Theoretical_Prob_dist);
sprintf('Lower bound value=%f',Lower_bound)
Thereotical_25th_percentile=chi2inv(0.25,4);
Thereotical_50th_percentile=chi2inv(0.50,4);
Thereotical_90th_percentile=chi2inv(0.90,4);
Estimated_25th_percentile=prctile(x,25);
Estimated_50th_percentile=prctile(x,50);
Estimated_90th_percentile=prctile(x,90);
sprintf('estimated 25th percentile using empirical distribution is=%f and Thereotical_25th percentile is=%f',Estimated_25th_percentile,Thereotical_25th_percentile)
sprintf('estimated 50th percentile using empirical distribution is=%f and Thereotical_50th percentile is=%f',Estimated_50th_percentile,Thereotical_50th_percentile)
sprintf('estimated 90th percentile using empirical distribution is=%f and Thereotical_90th percentile is=%f',Estimated_90th_percentile,Thereotical_90th_percentile)
end