function Proj3_Q3()
clc
Outcome_100=Generate_Random_Variable(100);
Outcome_1000=Generate_Random_Variable(1000);
Outcome_10000=Generate_Random_Variable(10000);

%below code plot Hist for 100 times
subplot(3,2,1)
hist(Outcome_100)
ylabel('frequency')
title('Histogram using 100 samples for N')
xlabel('Fewest # of samples required to get sum >4')
mean_100=mean(Outcome_100);
sprintf('mean of 100 samples for N=%f',mean_100)

%below code plot Hist for 1000 times
subplot(3,1,2)
hist(Outcome_1000)
ylabel('frequency')
title('Histogram using 1000 samples for N')
xlabel('Fewest # of samples required to get sum >4')
mean_1000=mean(Outcome_1000);
sprintf('mean of 1000 samples for N=%f',mean_1000)

%below code plot Hist for 10000 times
subplot(3,1,3)
hist(Outcome_10000)
ylabel('frequency')
title('Histogram using 10000 samples for N')
xlabel('Fewest # of samples required to get sum >4')
mean_10000=mean(Outcome_10000);
sprintf('mean of 10000 samples for N=%f',mean_10000)
end


function [N]=Generate_Random_Variable(No_of_times)
    sum_result=4;
    cal_sum_value=0;
    N=zeros(1,No_of_times);
    for i=1:No_of_times
       Length=1;
       Random_samples=[];
        while 1
            X=rand();
            Random_samples(Length)=X;
            
            if sum(Random_samples)>sum_result
                 N(i)=Length;
                 break;
            end
            Length=Length+1;
        end
   end
end