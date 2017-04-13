%------------------------------------------------------------------%
%Funtion name:count_No_of_Heads()
%Input parameters:None
%Output Parameters:None
%Defination:This function calls Flip_coins,Longest_run_heads(trial_outcomes)
%           and generates histogram of Bernoulli outcomes.
%--------------------------------------------------------------------%
function count_No_of_Heads()
 [No_of_Heads ,trial_outcomes]=Flip_coins();
 sprintf('Benoulli outcome(200 trials)=')
 disp(trial_outcomes)
 sprintf('NO of Heads in tossing a biased coin 200 times is=%d',No_of_Heads)
 Max_value=Longest_run_heads(trial_outcomes);
 sprintf('Longest run of heads is=%d',Max_value)
 histogram(trial_outcomes,2)
 xlabel('Heads or tails');
 ylabel('Frequency of heads and tails in 200 coin toss');
 title('Bernoulli outcomes of tossing a biased coin 200 times with P[HEADS]=0.8');
end
%------------------------------------------------------------------%
%Funtion name:Flip_coins()
%Input parameters:None
%Output Parameters:No_of_Heads ,trial_outcomes
%Defination:This function simulates tossing  of a biased coin 200 times,
%           count no. heads.         
%--------------------------------------------------------------------%

function [No_of_Heads ,trial_outcomes]=Flip_coins()
NO_of_Trials=200;
No_of_Heads=0;
trial_outcomes=zeros(1,NO_of_Trials);
outcomes=zeros(1,NO_of_Trials);
for iteration=1:NO_of_Trials
    % generate a number U[0,1] and threshold to fair Bernoulli trial
    outcomes(iteration) = rand();
    if outcomes(iteration) <0.800 %condition to test the outcome is head
        trial_outcomes(iteration)=1; %1 is considered as head and 0 is considered as tail
        No_of_Heads=No_of_Heads+1;
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
    Run_length(k)=Run_value;
end
%disp(Run_length);
Max_value=max(Run_length);
%disp(Max_value);
end 
