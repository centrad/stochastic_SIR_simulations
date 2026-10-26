## Work in progress

# lambda: vector probability of infecting an individual
# C contact rate in 1 time unit (number)
markov_sim_infectiousness <-function(T, lambda = c(), C, gamma, S0 = c(), I0 = c(), N) {
  t <- 0
  S <- S0
  I <- I0
  # Possible events
  eventLevels <- c("Ss->Is", "Sn->In", "Is->Rs", "In->Rn")
  # Initialize result
  events <- data.frame(t=0, 
                       Ss=S0[1],
                       Is=I0[1],
                       Sn=S0[2],
                       In=I0[2],
                       event=NA
                       )
  # Loop until we are past time T
  while (t < T & (as.integer(I[1] + I[2])>0)) {
    # Draw the waiting type for each possible event
    wait <-  rexp(4,c("Ss->Is" = (C*lambda[1]*S[1]*I[1]/N + 
                                  C*lambda[2]*S[1]*I[2]/N
                                ),
                      "Sn->In"=(C*lambda[2]*S[2]*I[2]/N +
                                C*lambda[1]*S[2]*I[1]/N
                                ),
                      "Is->Rs"=gamma*I[1],
                      "In->Rn"=gamma*I[2]
                      )
                  )
    # Which event occurs first
    i <- which.min(wait)
    # Advance Time
    t <- t+wait[i]
    # Update population according to the event type
    if (eventLevels[i] == "Ss->Is" & S[1] != 0) { 
      S[1] <- max(S[1]-1, 0) ; 
      I[1] <- I[1]+1
      }
    if (eventLevels[i] == "Sn->In" & S[2] != 0) { 
      S[2] <- max(S[2]-1, 0) ;
      I[2] <- I[2]+1
      }
    if (eventLevels[i] == "Is->Rs") { 
      I[1] <- max(I[1]-1, 0) 
      }
    if (eventLevels[i] == "In->Rn") { 
      I[2] <- max(I[2]-1, 0) 
      }
    #Store result
    events <- rbind(events,
                    c(t, as.integer(S[1]), as.integer(I[1]), as.integer(S[2]), as.integer(I[2]),i)
                    )
  }
  #Recode event type and return
  events$event <- factor(eventLevels[events$event], levels=eventLevels)
  return(events)
}




# lambda: probability of infecting an individual
# C contact matrix: contact rates between groups

markov_sim_contact <- function(T, lambda, C, gamma, S0 = c(), I0 = c(), N = c()) {
  t <- 0
  S <- S0
  I <- I0
  #Possible events
  eventLevels <- c("Ss->Is", "Sn->In", "Is->Rs", "In->Rn")
  #Initialize result
  events <- data.frame(t=0, 
                       Ss=S0[1],
                       Is=I0[1],
                       Sn=S0[2],
                       In=I0[2],
                       event=NA
                       )
  #Loop until we are past time T
  while (t < T & (as.integer(I[1] + I[2]) > 0)) {
    #Draw the waiting type for each possible event
    wait <-  rexp(4,
                  c("Ss->Is"=lambda*(C[1,1]*S[1]*I[1]/N[1] + 
                                     C[2,1]*S[1]*I[2]/N[1]),
                    "Sn->In"=lambda*(C[2,2]*S[2]*I[2]/N[2] + 
                                     C[1,2]*S[2]*I[1]/N[2]),
                    "Is->Rs"=gamma*I[1],
                    "In->Rn"=gamma*I[2])
                  )

    #Which event occurs first
    i <- which.min(wait)
    #Advance Time
    t <- t+wait[i]
    #Update population according to the event type
    if (eventLevels[i]=="Ss->Is" & S[1] != 0) { 
      S[1] <- max(S[1]-1,0) ; 
      I[1] <- I[1]+1
      }
    if (eventLevels[i]=="Sn->In" & S[2] != 0) { 
      S[2] <- max(S[2]-1,0) ; 
      I[2] <- I[2]+1
      }
    if (eventLevels[i]=="Is->Rs") { 
      I[1] <- max(I[1]-1,0) 
      }
    if (eventLevels[i]=="In->Rn") { 
      I[2] <- max(I[2]-1,0) 
      }
    #Store result
    events <- rbind(events,
                    c(t,as.integer(S[1]),as.integer(I[1]),as.integer(S[2]),as.integer(I[2]),i)
                    )
  }
  #Recode event type and return
  events$event <- factor(eventLevels[events$event], levels=eventLevels)
  return(events)
}




# table: dataframe with columns "contacts" (0,1,2,...) and "n" (number of individuals with contact number "contacts")
# fraction: a fraction of individuals with **lowest** contact number (0.9, 0.8,...,0.5) i.e. the socially non-superspreaders

# returns: the average contacts for the high-contacts group and for the low-contacts group

avg_contacts <- function(table, fraction){
  names(table) <- c("contacts", "n")
  individuals <- table %>% pull(n) %>% sum()
  n_perc <- fraction*individuals
  
  for(i in min(table$contacts):max(table$contacts)){
    tmp <- table %>% filter(contacts <= i) %>% pull(n) %>% sum()
    tmp1 <- table %>% filter(contacts <= i+1) %>% pull(n) %>% sum()
    if(tmp <n_perc & tmp1 >= n_perc){
      boundary <- i
      #print(boundary)
    }
  }
  
  # average contacts per individuals (in the n_perc group, i.e. socially non super-spreaders)
  average_low_den <- table %>% filter(contacts <= boundary) %>% pull(n) %>% sum()
  average_low_num <- table %>% filter(contacts <= boundary) %>% mutate(avg_num = contacts * n) %>% pull(avg_num) %>% sum()
  
  average_low <- average_low_num/average_low_den
  
  
  # average contacts per individuals (in the socially super-spreaders group)
  average_high_den <- table %>% filter(contacts > boundary) %>% pull(n) %>% sum()
  average_high_num <- table %>% filter(contacts > boundary) %>% mutate(avg_num = contacts * n) %>% pull(avg_num) %>% sum()
  
  average_high <- average_high_num/average_high_den
  
  return <- list("average_high_contact_group" = average_high, "average_low_contact_group" = average_low)
} 




# average_high: average contacts per individual belonging to the socially super-spreading group
# average_low: average contacts per individual belonging to the socially non super-spreading group
# type: specifies the type of social contacts
# fraction: the fraction of socially non super-spreaders in the population

# returns: a contact matrix C such that C_ij is the average number of contacts that an individual in group i has with an individual in group j. 
# in this framework: i=1 is the socially super-spreaders group.

contact_matrix <- function(average_high, average_low, type = c("assortative", "disassortative", "proportional"), fraction = NA){
  if(type == "assortative"){
    C <- diag(c(average_high, average_low))
  }
  if(type == "disassortative"){
    const <- fraction/(1-fraction)
    low_low <- max(0, average_low - average_high/const )
    off_diag <- average_low - low_low
     C <- t(matrix(c(average_high - off_diag * const, off_diag * const, off_diag, low_low), nrow=2))
  }
  if(type == "proportional"){
    p_high <- (1-fraction) * average_high/((1-fraction) * average_high+(fraction) * average_low)
    p_low <- fraction* average_low/((1-fraction) * average_high+(fraction) * average_low)
    C <- t(matrix(c(average_high *p_high, average_high * p_low, average_low * p_high, average_low * p_low), nrow=2))
  }
  return(C)
}







# params
# T max time
# lambda: vector of probability of infecting in an event of size c. length(lambda) = length(event_size). lambda[i] is the probability of infection by the infector participating in the event of size event_size[i] 
# C_mixing rate of events divided by population number (if all events are sized 2, it would correspond to C/2 in the heterogeneity (in infectiousness and in contacts where C was the contact rate per individual).

# 
#
# EX: rate of events: 10 / day. Population: 20 people. C_mixing: 10/20 = 0.5/day

# x  x  x  x  x  x  x  x  x  x   (10 people)
# |  |  |  |  |  |  |  |  |  |   ( meet )
# x  x  x  x  x  x  x  x  x  x   (10 people)

# C = 1 (each individual on average meets 1 person) 
# C = 1 = 2 C_mixing

# gamma recovery time
# S0 initial susceptibles
# I0 initial infectious
# N population number
# event_size: vector containing the range of event sizes. One could alternatively generate the event size according to a certain distribution
# prob_event_size: vector. prob_event_size[i] is the probability that the event has size event_size[i]


# returns: a number: 1 for major outbreaks, 0 otherwise

sim_mixing_groups_branching_process <- function(T, lambda = c(), C_mixing, gamma, S0 = c(), I0 = c(), N, event_size = c(), prob_event_size = c()){
  prob_event_size_biased <- prob_event_size * event_size/sum(event_size * prob_event_size)
  major <- 0
  t <- 0
  S <- S0
  I <- I0
  #Possible events
  eventLevels <- c("event", "I->R")
  
  #Initialize result
  events <- data.frame(t=0, S=S0,I=I0,event=NA)
  
  while(t < T & (as.integer(I)>0)){  
    
    #Draw the waiting type for each possible event
    wait <-  rexp(2,c("event"=(C_mixing*(sum(event_size * prob_event_size))*I),"I->R"=gamma*I))
    
    #Which event occurs first
    i <- which.min(wait)
    #Advance Time
    t <- t+wait[i]
    #Update population according to the event type
    if (eventLevels[i]=="event" & S != 0) { 
      # pick event size, with probability that the event containing that specific infected is of size c
      index <- sample(c(1:length(event_size)), size = 1, prob = prob_event_size_biased)
      # number of participants to the event
      Ne <- event_size[index]
      # probability of infection for that event size
      probability_inf <- lambda[index]
      # number of infected during the event
      inf_number <- rbinom(n = 1, size = Ne-1, prob = probability_inf)
      I <- I + inf_number
      S <- S - inf_number
    }
    if (eventLevels[i]=="I->R") { I <- max(I-1,0) }
    #Store result
    events <- rbind(events,c(t,as.integer(S),as.integer(I),i, NA))
    if((S0-S) > 0.1 * S0){
      major <- 1
      break
    }
  }
  #Recode event type and return
  events$event <- factor(eventLevels[as.integer(events$event)], levels=eventLevels)
  return(major)
}



############################################################################################################################################

# N = (NsH, NnH, NsL, NnL)
# ex: (N*0.1*0.1, N*0.9*0.1, N*0.1*0.9, N*0.9*0.9)
# lambda = (lambda_s, lambda_n)
# S = (SsH, SnH, SsL, SnL)
markov_sim_infectiousness_contact <- function(T, lambda = c(), C, gamma, S0 = c(), I0 = c(), N = c()) {
  
  t <- 0
  S <- S0
  I <- I0
  #Possible events
  eventLevels <- c("SsH->IsH", "SnH->InH", "SsL->IsL","SnL->InL", "IsH->RsH", "InH->RnH", "IsL->RsL","InL->RnL")
  #Initialize result
  events <- data.frame(t=0, SsH=S0[1],IsH=I0[1],SnH=S0[2],InH=I0[2],SsL=S0[3],IsL=I0[3], SnL=S0[3],InL=I0[3],event=NA)
  #Loop until we are past time T
  while (t < T & (as.integer(I[1] + I[2])>0)) {
    #Draw the waiting type for each possible event
    
    wait <-  rexp(8,c("SsH->IsH"=C[1,1]*S[1]/N[1] * (lambda[1]*I[1] + lambda[2]*I[2]) + C[2,1]*S[1]/N[1]*(lambda[1]*I[3] + lambda[2]*I[4]) ,
                      "SnH->InH"=C[1,1]*S[2]/N[2] * (lambda[1]*I[1] + lambda[2]*I[2]) + C[2,1]*S[2]/N[2]*(lambda[1]*I[3] + lambda[2]*I[4]), 
                      "SsL->IsL"=C[1,2]*S[3]/N[3] * (lambda[1]*I[1] + lambda[2]*I[2]) + C[2,2]*S[3]/N[3]*(lambda[1]*I[3] + lambda[2]*I[4]) ,
                      "SnL->InL"=C[1,2]*S[4]/N[4] * (lambda[1]*I[1] + lambda[2]*I[2]) + C[2,2]*S[4]/N[4]*(lambda[1]*I[3] + lambda[2]*I[4]) ,
                      "IsH->RsH"=gamma*I[1],
                      "InH->RnH"=gamma*I[2],
                      "IsL->RsL"=gamma*I[3],
                      "InL->RnL"=gamma*I[4]))
    # as above, with matrix notation
    # wait <-  rexp(4,c("Ss->Is"=t(S/N*(t(lambda*C) %*% I ))[1],"Sn->In"=t(S/N*(t(lambda*C) %*% I ))[2],"Is->Rs"=gamma*I[1],"In->Rn"=gamma*I[2]))
    
    #Which event occurs first
    i <- which.min(wait)
    #Advance Time
    t <- t+wait[i]
    #Update population according to the event type
    if (eventLevels[i]=="SsH->IsH" & S[1] != 0) { S[1] <- max(S[1]-1,0) ; I[1] <- I[1]+1}
    if (eventLevels[i]=="SnH->InH" & S[2] != 0) { S[2] <- max(S[2]-1,0) ; I[2] <- I[2]+1}
    if (eventLevels[i]=="SsL->IsL" & S[3] != 0) { S[3] <- max(S[3]-1,0) ; I[3] <- I[3]+1}
    if (eventLevels[i]=="SnL->InL" & S[4] != 0) { S[4] <- max(S[4]-1,0) ; I[4] <- I[4]+1}
    if (eventLevels[i]=="IsH->RsH") { I[1] <- max(I[1]-1,0) }
    if (eventLevels[i]=="InH->RnH") { I[2] <- max(I[2]-1,0) }
    if (eventLevels[i]=="IsL->RsL") { I[3] <- max(I[3]-1,0) }
    if (eventLevels[i]=="InL->RnL") { I[4] <- max(I[4]-1,0) }
    #Store result
    events <- rbind(events,c(t,as.integer(S[1]),as.integer(I[1]),as.integer(S[2]),as.integer(I[2]),as.integer(S[3]),as.integer(I[3]),as.integer(S[4]),as.integer(I[4]),i))
  }
  #Recode event type and return
  events$event <- factor(eventLevels[events$event], levels=eventLevels)
  return(events)
}


##################################################################################################################################################



sim_mixing_groups_branching_process_infectiousness <- function(T, lambda_s = c(), lambda_n = c(), C_mixing, gamma, S0 = c(), I0 = c(), event_size = c(), prob_event_size = c()){
  prob_event_size_biased <- prob_event_size * event_size/sum(event_size * prob_event_size)
  major <- 0
  t <- 0
  S <- S0
  I <- I0
  #Possible events
  eventLevels <- c("event_s", "event_n", "Is->Rs", "In->Rn")
  #Initialize result
  events <- data.frame(t=0, Ss=S0[1],Is=I0[1],Sn=S0[2],In=I0[2],event=NA)
  #Loop until we are past time T
  while (t < T & (as.integer(I[1] + I[2])>0)) {
    #Draw the waiting type for each possible event
    wait <-  rexp(4,c("event_s"=(C_mixing*(sum(event_size * prob_event_size))*I[1]),"event_n"=(C_mixing*(sum(event_size * prob_event_size))*I[2]),"Is->Rs"=gamma*I[1],"In->Rn"=gamma*I[2]))
    
    #Which event occurs first
    i <- which.min(wait)
    #Advance Time
    t <- t+wait[i]
    #Update population according to the event type
    if (eventLevels[i]=="event_s" & (S[1]+S[2]) != 0) { 
      # pick event size, with probability that the event containing that specific infected is of size c
      index <- sample(c(1:length(event_size)), size = 1, prob = prob_event_size_biased)
      # number of participants to the event
      Ne <- event_size[index]
      # probability of infection for that event size
      probability_inf <- lambda_s[index]
      # fraction of susceptible super-spreaders
      fraction <- S[1]/(S[1]+S[2])
      # pick number of super-spreaders susceptibles attending the meeting
      Ns <- rbinom(n = 1, size = Ne-1, prob = fraction)
      # number of infected during the event
      inf_number_s <- rbinom(n = 1, size = Ns, prob = probability_inf)
      inf_number_n <- rbinom(n = 1, size = Ne-1-Ns, prob = probability_inf)
      inf_number <- c(inf_number_s ,inf_number_n)
      I <- I + inf_number
      S <- S - inf_number
    }
    if (eventLevels[i]=="event_n" & (S[1]+S[2]) != 0) { 
      # pick event size, with probability that the event containing that specific infected is of size c
      index <- sample(c(1:length(event_size)), size = 1, prob = prob_event_size_biased)
      # number of participants to the event
      Ne <- event_size[index]
      # probability of infection for that event size
      probability_inf <- lambda_n[index]
      # fraction of susceptible super-spreaders
      fraction <- S[1]/(S[1]+S[2])
      # pick number of super-spreaders susceptibles attending the meeting
      Ns <- rbinom(n = 1, size = Ne-1, prob = fraction)
      # number of infected during the event
      inf_number_s <- rbinom(n = 1, size = Ns, prob = probability_inf)
      inf_number_n <- rbinom(n = 1, size = Ne-1-Ns, prob = probability_inf)
      inf_number <- c(inf_number_s ,inf_number_n)
      I <- I + inf_number
      S <- S - inf_number
    }
    if (eventLevels[i]=="Is->Rs") { I[1] <- max(I[1]-1,0) }
    if (eventLevels[i]=="In->Rn") { I[2] <- max(I[2]-1,0) }
    #Store result
    events <- rbind(events,c(t,as.integer(S[1]),as.integer(I[1]),as.integer(S[2]),as.integer(I[2]),i, NA))
    if((S0[1]-S[1]+ S0[2]-S[2]) > 0.1 * (S0[1]+S0[2])){
      major <- 1
      break
    }
  }
  #Recode event type and return
  events$event <- factor(eventLevels[as.integer(events$event)], levels=eventLevels)
  return(major)
}

##############################################################################################################



sim_mixing_groups_full_process <- function(T, lambda = c(), C_mixing, gamma, S0 = c(), I0 = c(), N, event_size = c(), prob_event_size = c()){

  t <- 0
  S <- S0
  I <- I0
  #Possible events
  eventLevels <- c("event", "I->R")
  
  #Initialize result
  events <- data.frame(t=0, S=S0,I=I0,event=NA)
  
  while(t < T & (as.integer(I)>0)){  
    
    if(S != 0){
      #Draw the waiting type for each possible event
      wait <-  rexp(2,c("event"=(N * C_mixing),"I->R"=gamma*I))
      
      #Which event occurs first
      i <- which.min(wait)
      #Advance Time
      t <- t+wait[i]
      #Update population according to the event type
      if (eventLevels[i]=="event" & S != 0) {
        # pick event size, with probability prob_event_size
        index <- sample(c(1:length(event_size)), size = 1, prob = prob_event_size)
        # number of participants to the event
        Ne <- event_size[index]
        # how many infectious people in the event?
        # fraction of infectious individuals
        fraction <- I/N
        # pick number of super-spreaders susceptibles attending the meeting
        Ni <- rbinom(n = 1, size = Ne, prob = fraction)
        Ns <- Ne-Ni
        # probability of infection for that event size
        probability_inf <- lambda[index]
        # number of infected during the event generated by each infectious
        inf_number <- ifelse(Ni != 0, rbinom(n = Ni, size = Ne-Ni, prob = probability_inf), 0)
        # non overlapping infections
        if(inf_number!=0){
          non_over = length(unique(sample(c(1:(Ne-Ni)), sum(inf_number), replace = TRUE)))
        }else(non_over = 0)
        
        I <- I + non_over
        S <- max(S - non_over,0)
       #print(I)
      }
      if (eventLevels[i]=="I->R") { I <- max(I-1,0) }
    } else if(S == 0){
      wait <- rexp(1, c("I->R"=gamma*I))
      t <- t+ wait
      I <- max(I-1, 0)
      i <- 2
      #print(I)
    }
   
    #Store result
    events <- rbind(events,c(t,as.integer(S),as.integer(I),i, NA))
  }
  #Recode event type and return
  events$event <- factor(eventLevels[as.integer(events$event)], levels=eventLevels)
  return(events)
}

