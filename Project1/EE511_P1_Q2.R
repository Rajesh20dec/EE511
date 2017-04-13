Flip_coins<-function(){
  NO_of_Heads<-0
  trial_outcomes<-sample(c(0:1),size=200,replace=TRUE,prob = c(0.2,0.8))
  cat("Outcome of 200 flips","\n", trial_outcomes,sep=" ")
  NO_of_Heads<-sum(trial_outcomes[which(trial_outcomes%%2==1)])
  Longest_run_heads(trial_outcomes)
  hist(trial_outcomes,breaks = 2,col="lightblue",
       main="Bernoulli outcomes of tossing a biased coin 200 times with P[HEADS]=0.8",
       xlab="Heads or tails",ylab="Frequency of heads and tails in 200 coin toss")
}


Longest_run_heads<-function(trial_outcomes){
  Max_value<-0
  Run_value<-0
  Run_length<-vector()
  k<-0
  for (i in c(1:length(trial_outcomes))){
    if (trial_outcomes[i]==1){
      Run_value<-Run_value+1
    }else{
      if ( Run_value>0){
        k<-k+1
        Run_length[k]<-Run_value
        Run_value<-0
      }
    }
  }
  if (Run_value >0){
    Run_length[k+1]<Run_value}
  Max_value<-max(Run_length)
  return(Max_value)
}