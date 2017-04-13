Project1_Q1<-function (){
  Repeat_1time()
  Repeat_20times()
  Repeat_100times()
  Repeat_200times()
  Repeat_1000times()
}


Repeat_1time<-function(){
  Return_value<-0
  Return_value<-Flip_coins()
  Max_value<-Longest_run_heads(Return_value$trial_outcomes)
  colors<-c("red","green")
  hist(Return_value$trial_outcomes,breaks = 2,col=colors,
      main="Tossing of fair coin 50 times",
      xlab="outcome",ylab="frequency of outcome")
  
}


Repeat_20times<-function(){
  Freq_of_Heads=rep(0,20)
  Return_value<-0
  for (i in c(1:20)){
  Return_value<-Flip_coins()
  Freq_of_Heads[i]=Return_value$NO_of_Heads
  }

  cat("Freq_of_Heads","\n", Freq_of_Heads,sep=" ","\n")
  hist(Freq_of_Heads,col="lightblue",
       main="Tossing of fair coin 20*50 times",
       xlab="No. of Heads in 50 flips",ylab="Frequency of nos. of heads in 50 flips")
}


Repeat_100times<-function(){
  Freq_of_Heads=rep(0,100)
  Return_value<-0
  for (i in c(1:100)){
    Return_value<-Flip_coins()
    Freq_of_Heads[i]=Return_value$NO_of_Heads
  }
  cat("Freq_of_Heads","\n", Freq_of_Heads,sep=" ","\n")
  hist(Freq_of_Heads,col="lightblue",
       main="Tossing of fair coin 100*50 times",
       xlab="No. of Heads in 50 flips",ylab="Frequency of nos. of heads in 50 flips")
 
}



Repeat_200times<-function(){
  Freq_of_Heads=rep(0,200)
  Return_value<-0
  for (i in c(1:200)){
    Return_value<-Flip_coins()
    Freq_of_Heads[i]=Return_value$NO_of_Heads
  }
  cat("Freq_of_Heads","\n", Freq_of_Heads,sep=" ","\n")
  hist(Freq_of_Heads,col= "lightblue",
       main="Tossing of fair coin 200*50 times",
       xlab="No. of Heads in 50 flips",ylab="Frequency of nos. of heads in 50 flips")
  #hist(Freq_of_Heads,col="blue")
}


Repeat_1000times<-function(){
  Freq_of_Heads<-rep(0,1000)
  Return_value<-0
  for (i in c(1:1000)){
    Return_value<-Flip_coins()
    Freq_of_Heads[i]=Return_value$NO_of_Heads
  }
  cat("Freq_of_Heads","\n", Freq_of_Heads,sep=" ","\n")
  hist(Freq_of_Heads,col= "lightblue",
       main="Tossing of fair coin 1000*50 times",
       xlab="No. of Heads in 50 flips",ylab="Frequency of nos. of heads in 50 flips")
 # hist(Freq_of_Heads,col="blue")
}


Flip_coins<-function(){
  NO_of_Heads<-0
  trial_outcomes<-sample(c(0,1),50,replace=TRUE)
  #print(trial_outcomes)
  NO_of_Heads<-sum(trial_outcomes[which(trial_outcomes%%2==1)])
  My_list<-list("trial_outcomes"=trial_outcomes,"NO_of_Heads"=NO_of_Heads)
  return(My_list)
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
 


