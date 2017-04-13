%------------------------------------------------------------------%
%Funtion name:Project1_Q1()
%Input parameters:None
%Output Parameters:None
%Defination:Each problem of question 1 is implemented in separte functions
%            and this function calls these functions.
%--------------------------------------------------------------------%
function Project1_Q1()
Repeat_1time();
Repeat_20times();
Repeat_100times();
Repeat_200times();
Repeat_1000times();
end
%------------------------------------------------------------------%
%Funtion name:Repeat_1time()
%Input parameters:None
%Output Parameters:None
%Defination:This functions call Flip_coin()
%           and Longest_run_heads(trial_outcomes)
%           and generates Histogram of Bernoull outcomes.
%--------------------------------------------------------------------%
function Repeat_1time()
[No_of_Heads ,trial_outcomes]=Flip_coins();
%disp(No_of_Heads)
%disp(trial_outcomes)
sprintf('No.of heads in trial %d',No_of_Heads)
Max_value=Longest_run_heads(trial_outcomes);
sprintf('Benoulli outcome(50 trials)=')
disp(trial_outcomes)
sprintf('Longest runs of heads is= %d',Max_value)
figure(1);
histogram(trial_outcomes,2)
xlabel('outcome');
ylabel('frequency of outcome')
title('Tossing of fair coin 50 times');
end
%------------------------------------------------------------------%
%Funtion name:Repeat_20times()
%Input parameters:None
%Output Parameters:None
%Defination:This functions call Flip_coin() 20 times(1 time cosnists of 50
%           toss of coins) and calculate nos. of heads in each experiment
%           of 50 trials of toss and generates Histogram of nos of heads 
%            in 50 filps.
%--------------------------------------------------------------------%

function Repeat_20times()
Freq_of_Heads=zeros(1,20);
for i=1:20
    [No_of_Heads ,trial_outcomes]=Flip_coins();
    Freq_of_Heads(i)=No_of_Heads;
    
end
%disp(sort(Freq_of_Heads))
figure(2);
histogram(Freq_of_Heads)
xlabel('No. of Heads in 50 flips');
ylabel('Frequency of nos. of heads in 50 flips')
title('Tossing of fair coin 20*50 times');
end
%------------------------------------------------------------------%
%Funtion name:Repeat_100times()
%Input parameters:None
%Output Parameters:None
%Defination:This functions call Flip_coin() 100 times(1 time cosnists of 50
%           toss of coins) and calculate nos. of heads in each experiment
%           of 50 trials of toss and generates Histogram of nos of heads 
%            in 50 filps.
%--------------------------------------------------------------------%

function Repeat_100times()
Freq_of_Heads=zeros(1,100);
for i=1:100
    [No_of_Heads ,trial_outcomes]=Flip_coins();
    Freq_of_Heads(i)=No_of_Heads;
end
%disp(sort(Freq_of_Heads))
figure(3);
histogram(Freq_of_Heads)
xlabel('No. of Heads in 50 flips');
ylabel('Frequency of nos. of heads in 50 flips')
title('Tossing of fair coin 100*50 times');
end
%------------------------------------------------------------------%
%Funtion name:Repeat_200times()
%Input parameters:None
%Output Parameters:None
%Defination:This functions call Flip_coin() 200 times(1 time cosnists of 50
%           toss of coins) and calculate nos. of heads in each experiment
%           of 50 trials of toss and generates Histogram of nos of heads 
%            in 50 filps.
%--------------------------------------------------------------------%
function Repeat_200times()
Freq_of_Heads=zeros(1,200);
for i=1:200
    [No_of_Heads ,trial_outcomes]=Flip_coins();
    Freq_of_Heads(i)=No_of_Heads;
end
%disp(sort(Freq_of_Heads))
figure(4);
histogram(Freq_of_Heads)
xlabel('No. of Heads in 50 flips');
ylabel('Frequency of nos. of heads in 50 flips')
title('Tossing of fair coin 200*50 times');
end
%------------------------------------------------------------------%
%Funtion name:Repeat_1000times()
%Input parameters:None
%Output Parameters:None
%Defination:This functions call Flip_coin() 1000 times(1 time cosnists of 50
%           toss of coins) and calculate nos. of heads in each experiment
%           of 50 trials of toss and generates Histogram of nos of heads 
%            in 50 filps.
%--------------------------------------------------------------------%
function Repeat_1000times()
Freq_of_Heads=zeros(1,1000);
for i=1:1000
    [No_of_Heads ,trial_outcomes]=Flip_coins();
    Freq_of_Heads(i)=No_of_Heads;
end

%disp(sort(Freq_of_Heads))
figure(5);
histogram(Freq_of_Heads)
xlabel('No. of Heads in 50 flips');
ylabel('Frequency of nos. of heads in 50 flips')
title('Tossing of fair coin 1000*50 times');
end
%------------------------------------------------------------------%
%Funtion name:Flip_coins()
%Input parameters:None
%Output Parameters:No_of_Heads ,trial_outcomes
%Defination:This function simulates tossing  of a coin 50 times,
%           count no. heads.
%           
%--------------------------------------------------------------------%
    
function [No_of_Heads ,trial_outcomes]=Flip_coins()
NO_of_Trials=50;
No_of_Heads=0;
trial_outcomes=zeros(1,NO_of_Trials);
outcomes=zeros(1,NO_of_Trials);
for iteration=1:NO_of_Trials
    % generate a number U[0,1] and threshold to fair Bernoulli trial
        outcomes(iteration) = rand();
    if  outcomes(iteration) >0.5000  %condition to test the outcome is head
        trial_outcomes(iteration)=1; %1 is considred as head and 0 is considered as tail
        No_of_Heads=No_of_Heads+1;    % sum up the heads count
    end
end
end
%------------------------------------------------------------------%
%Funtion name:Longest_run_heads(trial_outcomes)
%Input parameters:trial_outcomes
%Output Parameters:Max_value.
%Defination:This function takes bernoulli outcomes as input and compute
%           longest run of heads
%--------------------------------------------------------------------%

function Max_value=Longest_run_heads(trial_outcomes)
Run_length=[];
k=0;
Run_value=0;
Max_value=0;
for i=1:length(trial_outcomes)
    if trial_outcomes(i)==1
        Run_value=Run_value+1;
    else
        if Run_value>0
            k=k+1;
            Run_length(k)=Run_value;
            Run_value=0;
        end
        
    end
end
if Run_value>0
    Run_length(k+1)=Run_value;
end
%disp(Run_length);
Max_value=max(Run_length);
%disp(Max_value);
end  

     
        
   