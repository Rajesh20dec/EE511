function Proj_P4_Q3()
clc
clear all
Waiting_time=[79,54,74,62,85,55,88,85,51,85,54,84,78,47,83];
Mean_15=mean(Waiting_time);
size_of_sample=15;
No_of_samples=1000;
mean_bar=[];
for i=1:No_of_samples
    for j=1:size_of_sample
        %collecting samples of size =5 from  waiting times range taken as input as 
        sample(j)=randi([44 125]);
    end
mean_of_collected_samples(i)=mean(sample);
STD_of_collected_samples(i)=std(sample);
STD_error(i)= STD_of_collected_samples(i)./(sqrt(size_of_sample));
T_distribution(i)=(mean_of_collected_samples(i)-Mean_15)./STD_error(i);
Margin_of_error(i)=STD_error(i)*T_distribution(i);
end
Mean_of_margin_error= mean(Margin_of_error);
CI__Upper_value=Mean_15+Mean_of_margin_error;
CI_Lower_value=Mean_15-Mean_of_margin_error;
%below line will calculate theroretical Confidence interval 
Theroretical_CI_Range=bootci(No_of_samples,@mean,mean_of_collected_samples);
sprintf('Extimated  value of confidence Interval is[%f %f]',CI_Lower_value,CI__Upper_value)
sprintf('Theoretical value of confidence interval is [%f %f]',Theroretical_CI_Range(1),Theroretical_CI_Range(2))
end