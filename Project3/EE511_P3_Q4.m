function Proj_P3_Q4()
clc
Generate_Xk();
Generate_min_N60()
end

function Generate_Xk()
p=1/sum(1./[1:60]);%this calculates p
sprintf('normalisation value of p =%f',p)
Prob=p./[1:1:60];  %this calculates probality distribution
Samples=rand(1,1000);
random_sample=[];
Prob_i_minus_1=0;
prob_i=0;
for i=1:60
    counter=0;
    prob_i=prob_i+Prob(i);
    for j=1:1000
        
        if Samples(j)<prob_i && Samples(j)>=Prob_i_minus_1 %this checks condition  summation of pi i=0 through i-1 ?U < summation of pi i=0 to i=j
            counter=counter+1;
        end   
    end
    random_sample(i)=counter;
    Prob_i_minus_1=prob_i;
end

figure(1)
bar(1:60,random_sample)
title('histogram of Xk')
xlabel('Xk=i')
ylabel('Frequency of Xk=i')
end

function Generate_min_N60()
p=1/sum(1./[1:60]);%this calculates p
Prob_N60=p/60;    %Prob. of P60
total_experiment=10000;
Min_Sample=[];

for i=1:total_experiment
    iteration=0;
    while rand>=Prob_N60
            iteration=iteration+1;
    end
    Min_Sample(i)=iteration;
end
figure(2)
    histogram(Min_Sample);
    title('Histogram of min # of samples require to get N60  with # of experiment =10000') 
    xlabel('Min. # of samples to get N60')
    ylabel('Frequency')
    sample_mean=mean(Min_Sample);
    sample_variance=var(Min_Sample);
    theoretical_mean=1/Prob_N60;
    theoretical_variance=(1-Prob_N60)/(Prob_N60*Prob_N60);
    sprintf('sample mean =%f and theoretical_mean=%f',sample_mean,theoretical_mean)
    sprintf('sample variance =%f and theoretical_variance=%f',sample_variance,theoretical_variance)
end
