function simulate_1hr_arrival()
clc;
lamda=120;
No_of_experiment=1000;
No_of_car_arrival_Bernoulli=zeros(1,No_of_experiment);
No_of_car_arrival_Inverse=zeros(1,No_of_experiment);


%Below code generates ramdom samples using subinterval method

NO_of_subinterval=100000;
p=lamda/NO_of_subinterval;
    for i =1:No_of_experiment
        x=rand(NO_of_subinterval,1);
        bernoulliTrials=x<p;
        No_of_car_arrival_Bernoulli(i) = sum(bernoulliTrials);
    end

    
%Below code genearates Ramdom samples using Inverse transorm method

    for i=1 :No_of_experiment
        flag=-1;
        j=0;
        U=rand();
        P=0;
        while flag~=0
            P=poisspdf(j,120)+P;
            if U<P
                flag=0;
            end
            j=j+1;
        end
        No_of_car_arrival_Inverse(i)=j-1;
    end

%below code generates Possion distriobution 1000 times
    for i=1:1000
        outcomes(i)=poisspdf(i,120);
    end
theoretical=1:1:1000;

%below code is used to plot hist of bernoulli overlays 
%with therotical p.m.f

figure(1)
yyaxis left
h=stem(theoretical,outcomes);
h.Color = 'red';
yyaxis right
hist(No_of_car_arrival_Inverse);
title('Plot of Random variable using Inverse transform and theoretical p.m.f')
xlabel('# of experiment')
yyaxis left
ylabel('P.m.f value')
yyaxis right
ylabel('Frequncy of # of car arrival per hour')
legend({'p.m.f','Using Inverse transform'},'FontSize',20)

%below code is used to plot hist of inverse transform overlays 
%with therotical p.m.f

figure(2)
yyaxis left
h=stem(theoretical,outcomes);
h.Color = 'red';
yyaxis right
hist(No_of_car_arrival_Bernoulli);
title('Plot of Random variable using Bernoulli trial and theoretical p.m.f')
xlabel('# of experiment')
yyaxis left
ylabel('P.m.f value')
yyaxis right
ylabel('Frequncy of # of car arrival per hour')
legend({'p.m.f','Using Bernoulli'},'FontSize',20)

end
