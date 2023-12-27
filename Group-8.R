# Importing package 
library(deSolve) # solve ode 
library(ggplot2) # plot graph 
library(tidyr) # plot graph

#########################################################
### FUNCTION TO VISUALIZE THE MODEL 
plot_result <- function(model){
  model_df <- as.data.frame(model)
  # Convert the dataframe to long format for ggplot
  model_df_long <- gather(model_df, key = "Group", value = "Value", -time)
  
  # Set the order of levels for the "Group" factor
  model_df_long$Group <- factor(model_df_long$Group, levels = c("S", "E", "I", "R", "D"))
  
  # Plotting using ggplot
  ggplot(model_df_long, aes(x = time, y = Value, color = Group)) +
    geom_line() +
    labs(title = "Evolution of epidemic",
         x = "Time",
         y = "Population") +
    scale_color_manual(values = c("blue", "yellow", "red", "green", "black")) +
    theme_minimal()
}
###
#########################################################
## SEIR

SEIR <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+E+I+R
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I
    dR <- gamma*I
    
    return(list(c(dS, dE, dI, dR)))
  })
}

## parameters 
params <- c(
  beta=0.5, 
  sigma=0.25, 
  gamma=0.2
)
###

### initial state
initial_state <- c(S=999997, E=3, I=0, R=0)
times <- 0:365
### 


model <- ode(initial_state, times, SEIR, params)

print(summary(model))

plot_result(model)


######################################################

#### Introduing new compartment D(Deceased)
# without intervention 

SEIRD <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+E+I+R+D
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I
    dD <- mu*I
    
    return(list(c(dS, dE, dI, dR, dD)))
  })
}

# parameters
params <- c(beta=0.5, 
            sigma=0.25, 
            gamma=0.2, 
            mu=0.001)

# Initial State
initial_state <- c(S=999997, E=3, I=0, R=0, D=0)
times <- 0:365


model_2 <- ode(initial_state, times, SEIRD, params)
summary(model_2)
plot_result(model_2)

######################################################

### new variant

# Initial State
initial_state <- c(S=999997, E=3, I=0, R=0, D=0)
times <- 0:365

params <- c(beta=0.5, 
            sigma=0.25, 
            gamma=0.2, 
            mu=0.001,
            start = 90)

SEIRD_new_variant <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    gamma = ifelse(
      (time <= start),
      0.2, 0.07
    )
    
    mu = ifelse(
      (time <= start),
      0.001, 0.015
    )
    
    N <- S+E+I+R+D
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I
    dD <- mu*I
    
    return(list(c(dS, dE, dI, dR, dD)))
  })
}

model_3 <- ode(initial_state, times, SEIRD_new_variant, params)
summary(model_3)
plot_result(model_3)

###################################################################

### with intervention
SEIRD_lockdown <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    beta = ifelse(
      (time <= start_lockdown || time >= end_lockdown),
      0.5, 0.1
    )
    
    N <- S+E+I+R+D
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I
    dD <- mu*I
    
    return(list(c(dS, dE, dI, dR, dD)))
  })
}

params <- c(
  sigma=0.25,
  gamma=0.2,
  mu=0.001,
  start_lockdown=90,
  end_lockdown=150
)

initial_state <- c(S=999997, E=3, I=0, R=0, D=0)

times <- 0:365

model_4 <- ode(initial_state, times, SEIRD_lockdown, params)
summary(model_4)
plot_result(model_4)

#########################################
# vaccination

SEIRD_vaccine <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+E+I+R+D
    dS <- -(beta*S*I)/N - v*S/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - (gamma+mu)*I
    dR <- gamma*I + v*S/N
    dD <- mu*I
    
    return(list(c(dS, dE, dI, dR, dD)))
  })
}

# parameters
params <- c(beta=0.5, 
            sigma=0.25, 
            gamma=0.2, 
            mu=0.001,
            v = 3000)

# Initial State
initial_state <- c(S=999997, E=3, I=0, R=0, D=0)
times <- 0:365

model_5 <- ode(initial_state, times, SEIRD_vaccine, params)
summary(model_5)
plot_result(model_5)
  
##################################
