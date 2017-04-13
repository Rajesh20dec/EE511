Flip_coins<-function(){
  Run_length<-vector()
  trial_outcomes<-sample(c(0:1),size=100,replace=TRUE,prob = c(0.5,0.5))
  Run_length<-Heads_run_lengths(trial_outcomes)
  cat("Outcome of 100 flips","\n", trial_outcomes,sep=" ")
  cat("Run Lenths of Heads=","\n", Run_length,sep=" ")
  #print(Run_length)
  hist(Run_length,col="green")
  hist(Run_length,col="lightblue",
       main="Histogram showing the heads run length",
       xlab="Heads or tails no.",ylab="Frequency")

  
}



Heads_run_lengths<-function(trial_outcomes){
  
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
    return(Run_length)
}