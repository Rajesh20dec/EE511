function Sample_Uniformly()
NO.of_Times=input('Enter how many times sampling is required:');
Lower_Lim=-3;
Upper_Lim=2;
Outcomes=Lower_Lim+(Upper_Lim-Lower_Lim)*rand(1,NO.of_Times);%a+(b-a)*rand()
%sprintf('Outcome of experiment=')
%disp(Outcomes)
xlim([-3 2]);
hist(Outcomes)
xlabel('outcomes');
ylabel('Frequency');
title('uniform sample');
Sample_Mean=mean(Outcomes);
Therotical_Mean=(Lower_Lim+Upper_Lim)/2;
sprintf('Sample Mean =%f and Therotical_Mean=%f',Sample_Mean,Therotical_Mean)
Sample_Variance=var(Outcomes);
Therotical_Variance=(Upper_Lim-Lower_Lim)^2/12;
sprintf('Sample Variance = %f and Therotical_Variance=%f',Sample_Variance,Therotical_Variance)
%Bootstrap_samples=bootstrp(1000,@mean,Outcomes);
Confidence_Interval_mean=bootci(1000,@mean,Outcomes);
sprintf('Bootstrap confidenct interval for Sample mean=')
disp(Confidence_Interval_mean)
Confidence_Interval_Variance=bootci(1000,@std,Outcomes);
sprintf('Bootstrap confidenct interval for Sample variance=%')
disp(Confidence_Interval_Variance)
end


