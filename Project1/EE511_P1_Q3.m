%------------------------------------------------------------------%
%Funtion name:show_the_heads_run_length()
%Input parameters:None
%Output Parameters:None
%Defination:This function shows head run  lenghts
%                   
%--------------------------------------------------------------------%
function show_the_heads_run_length()
trial_outcomes =Flip_coins();
Length=Heads_run_lengths(trial_outcomes);
histogram(Length)
xlabel('Heads or tails no.');
ylabel('Frequency');
title('Histogram showing the heads run length');
end
%------------------------------------------------------------------%
%Funtion name:Flip_coins()
%Input parameters:None
%Output Parameters:trial_outcomes
%Defination:his function simulates tossing  of a coin 100 times,
%           count no. heads.        
%--------------------------------------------------------------------%
function trial_outcomes =Flip_coins()
NO_of_Trials=100;
trial_outcomes=zeros(1,NO_of_Trials);
outcomes=zeros(1,NO_of_Trials);
for iteration=1:NO_of_Trials
    % generate a number U[0,1] and threshold to fair Bernoulli trial
    outcomes(iteration) = rand();
    if outcomes(iteration) <0.500 %condition to test the outcome is head 
        trial_outcomes(iteration)=1; %1 is considered head and 0 tail
    end
end
end
%------------------------------------------------------------------%
%Funtion name:Heads_run_lenths()
%Input parameters:trial_outcomes
%Output Parameters:Length
%Defination:This function calculate heads run length 
%                    
%--------------------------------------------------------------------%
function Length=Heads_run_lengths(trial_outcomes)
Length=[];
k=0;
Run_value=0;
for i=1:length(trial_outcomes)
    if (trial_outcomes(i)==1)
        Run_value=Run_value+1;
    else
        if Run_value>0
            k=k+1;
            Length(k)=Run_value;
            Run_value=0;
         end
   end
end
if Run_value>0
 Length(k+1)=Run_value;
end
sprintf('Benoulli outcome(100 trials)=')
disp(trial_outcomes)
%disp(Length)
%Length=Length(1:k-1);
sprintf('Heads length=')
disp(Length)
end