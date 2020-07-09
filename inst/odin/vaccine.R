################################################################################
### Vaccine model: deterministic ###############################################
################################################################################

### Initial setup ##############################################################
dt <- user() # Specified timestep
N_age <- user()
N_vaccine <- user()
time <- t
output(time) <- TRUE
################################################################################

### S: susceptible #############################################################
dim(S) <- c(N_age, N_vaccine)

S_0[,] <- user()
dim(S_0) <- c(N_age, N_vaccine)
initial(S[,]) <- S_0[i,j]

deriv(S[,1]) <- (gamma_R * R2[i,j]) - (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j]) - (vr * vaccination_target[i] * S[i,j])
deriv(S[,2]) <- (vr * vaccination_target[i] * S[i,j-1]) + (gamma_R * R2[i,j]) - (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j]) - (gamma_vaccine[j] * S[i,j])
deriv(S[,3:N_vaccine]) <- (gamma_vaccine[j-1] * S[i,j-1]) + (gamma_R * R2[i,j]) - (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j]) - (gamma_vaccine[j] * S[i,j])

output(Sout[]) <- sum(S[i,])
dim(Sout) <- N_age
################################################################################

### E (E1 & E2): Latent ########################################################
dim(E1) <- c(N_age, N_vaccine)
dim(E2) <- c(N_age, N_vaccine)

E1_0[,] <- user()
dim(E1_0) <- c(N_age, N_vaccine)
initial(E1[,]) <- E1_0[i,j]

E2_0[,] <- user()
dim(E2_0) <- c(N_age, N_vaccine)
initial(E2[,]) <- E2_0[i,j]

gamma_E <- user() # rate of progression through latent infection

deriv(E1[,1]) <- (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j]) - (gamma_E * E1[i,j]) - (vr * vaccination_target[i] * E1[i,j])
deriv(E1[,2]) <- (vr * vaccination_target[i] * E1[i,j-1]) + (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j]) - (gamma_E * E1[i,j]) - (gamma_vaccine[j] * E1[i,j])
deriv(E1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * E1[i,j-1]) + (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j]) - (gamma_E * E1[i,j]) - (gamma_vaccine[j] * E1[i,j])

deriv(E2[,1]) <- (gamma_E * E1[i,j]) - (gamma_E * E2[i,j]) - (vr * vaccination_target[i] * E2[i,j])
deriv(E2[,2]) <- (vr * vaccination_target[i] * E2[i,j-1]) + (gamma_E * E1[i,j]) - (gamma_E * E2[i,j]) - (gamma_vaccine[j] * E2[i,j])
deriv(E2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * E2[i,j-1]) + (gamma_E * E1[i,j]) - (gamma_E * E2[i,j]) - (gamma_vaccine[j] * E2[i,j])

output(Eout[]) <- sum(E1[i,]) + sum(E2[i,])
dim(Eout) <- N_age
################################################################################

### IMild: Unhospitalised infection ############################################
dim(IMild) <- c(N_age, N_vaccine)

IMild_0[,] <- user()
dim(IMild_0) <- c(N_age, N_vaccine)
initial(IMild[,]) <- IMild_0[i,j]

gamma_IMild <- user() # rate of progression from mild infection to recovery

deriv(IMild[,1]) <- (gamma_E * E2[i,j] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j])
deriv(IMild[,2]) <- (gamma_E * E2[i,j] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * IMild[i,j])
deriv(IMild[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMild[i,j-1]) + (gamma_E * E2[i,j] * (1 - prob_hosp[i,j])) - (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * IMild[i,j])

output(IMildout[]) <- sum(IMild[i,])
dim(IMildout) <- N_age
################################################################################

### R: (R1 & R2): Recovered ####################################################
dim(R1) <- c(N_age, N_vaccine)
dim(R2) <- c(N_age, N_vaccine)

R1_0[,] <- user()
dim(R1_0) <- c(N_age, N_vaccine)
initial(R1[,]) <- R1_0[i,j]

R2_0[,] <- user()
dim(R2_0) <- c(N_age, N_vaccine)
initial(R2[,]) <- R2_0[i,j]

gamma_R <- user() # rate of progression through recovered compartment (loss of naturally acquired immunity)

deriv(R1[,1]) <- (gamma_rec * IRec2[i,j]) + (gamma_IMild * IMild[i,j]) + (gamma_get_ox_survive * IOxGetLive2[i,j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_R * R1[i,j]) - (vr * vaccination_target[i] * R1[i,j])
deriv(R1[,2]) <- (vr * vaccination_target[i] * R1[i,j-1]) + (gamma_rec * IRec2[i,j]) + (gamma_IMild * IMild[i,j]) + (gamma_get_ox_survive * IOxGetLive2[i,j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_R * R1[i,j]) - (gamma_vaccine[j] * R1[i,j])
deriv(R1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * R1[i,j-1]) + (gamma_rec * IRec2[i,j]) + (gamma_IMild * IMild[i,j]) + (gamma_get_ox_survive * IOxGetLive2[i,j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_R * R1[i,j]) - (gamma_vaccine[j] * R1[i,j])

deriv(R2[,1]) <- (gamma_R * R1[i,j]) - (gamma_R * R2[i,j]) - (vr * vaccination_target[i] * R2[i,j])
deriv(R2[,2]) <- (vr * vaccination_target[i] * R2[i,j-1]) + (gamma_R * R1[i,j]) - (gamma_R * R2[i,j]) - (gamma_vaccine[j] * R2[i,j])
deriv(R2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * R2[i,j-1]) + (gamma_R * R1[i,j]) - (gamma_R * R2[i,j]) - (gamma_vaccine[j] * R2[i,j])

output(Rout[]) <- sum(R1[i,]) + sum(R2[i,])
dim(Rout) <- N_age
################################################################################

### ICase (ICase1 & ICase2): To-be hospitalised infection ######################
dim(ICase1) <- c(N_age, N_vaccine)
dim(ICase2) <- c(N_age, N_vaccine)

ICase1_0[,] <- user()
dim(ICase1_0) <- c(N_age, N_vaccine)
initial(ICase1[,]) <- ICase1_0[i,j]

ICase2_0[,] <- user()
dim(ICase2_0) <- c(N_age, N_vaccine)
initial(ICase2[,]) <- ICase2_0[i,j]

gamma_ICase <- user() # rate of progression from symptom onset to requiring hospitalisation

deriv(ICase1[,1]) <- (gamma_E * E2[i,j] * prob_hosp[i,j]) - (gamma_ICase * ICase1[i,j])
deriv(ICase1[,2]) <- (gamma_E * E2[i,j] * prob_hosp[i,j]) - (gamma_ICase * ICase1[i,j]) - (gamma_vaccine[j] * ICase1[i,j])
deriv(ICase1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * ICase1[i,j-1]) + (gamma_E * E2[i,j] * prob_hosp[i,j]) - (gamma_ICase * ICase1[i,j]) - (gamma_vaccine[j] * ICase1[i,j])

deriv(ICase2[,1]) <- (gamma_ICase * ICase1[i,j]) - (gamma_ICase * ICase2[i,j])
deriv(ICase2[,2]) <- (gamma_ICase * ICase1[i,j]) - (gamma_ICase * ICase2[i,j]) - (gamma_vaccine[j] * ICase2[i,j])
deriv(ICase2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * ICase2[i,j-1]) + (gamma_ICase * ICase1[i,j]) - (gamma_ICase * ICase2[i,j]) - (gamma_vaccine[j] * ICase2[i,j])

output(ICaseout[]) <- sum(ICase1[i,]) + sum(ICase2[i,])
dim(ICaseout) <- N_age
################################################################################

### IOxGetLive (IOxGetLive1 & IOxGetLive2): Get oxygen, go on to survive #######
dim(IOxGetLive1) <- c(N_age, N_vaccine)
dim(IOxGetLive2) <- c(N_age, N_vaccine)

IOxGetLive1_0[,] <- user()
dim(IOxGetLive1_0) <- c(N_age, N_vaccine)
initial(IOxGetLive1[,]) <- IOxGetLive1_0[i,j]

IOxGetLive2_0[,] <- user()
dim(IOxGetLive2_0) <- c(N_age, N_vaccine)
initial(IOxGetLive2[,]) <- IOxGetLive2_0[i,j]

gamma_get_ox_survive <- user() # rate of progression through requiring oxygen compartment conditional on getting oxygen and surviving

deriv(IOxGetLive1[,1]) <- ((1 - prob_non_severe_death_treatment[i]) * number_get_Ox[i,j]) - (gamma_get_ox_survive * IOxGetLive1[i,j])
deriv(IOxGetLive1[,2]) <- ((1 - prob_non_severe_death_treatment[i]) * number_get_Ox[i,j]) - (gamma_get_ox_survive * IOxGetLive1[i,j]) - (gamma_vaccine[j] * IOxGetLive1[i,j])
deriv(IOxGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxGetLive1[i,j-1]) + ((1 - prob_non_severe_death_treatment[i]) * number_get_Ox[i,j]) - (gamma_get_ox_survive * IOxGetLive1[i,j]) - (gamma_vaccine[j] * IOxGetLive1[i,j])

deriv(IOxGetLive2[,1]) <- (gamma_get_ox_survive * IOxGetLive1[i,j]) - (gamma_get_ox_survive * IOxGetLive2[i,j])
deriv(IOxGetLive2[,2]) <- (gamma_get_ox_survive * IOxGetLive1[i,j]) - (gamma_get_ox_survive * IOxGetLive2[i,j]) - (gamma_vaccine[j] *  IOxGetLive2[i,j])
deriv(IOxGetLive2[,3:N_vaccine]) <- (gamma_vaccine[j-1] *  IOxGetLive2[i,j-1]) + (gamma_get_ox_survive * IOxGetLive1[i,j]) -  (gamma_get_ox_survive * IOxGetLive2[i,j]) - (gamma_vaccine[j] *  IOxGetLive2[i,j])
################################################################################

### IOxGetDie (IOxGetDie1 & IOxGetDie2): Get oxygen go on to die ###############
dim(IOxGetDie1) <- c(N_age, N_vaccine)
dim(IOxGetDie2) <- c(N_age, N_vaccine)

IOxGetDie1_0[,] <- user()
dim(IOxGetDie1_0) <- c(N_age, N_vaccine)
initial(IOxGetDie1[,]) <- IOxGetDie1_0[i,j]

IOxGetDie2_0[,] <- user()
dim(IOxGetDie2_0) <- c(N_age, N_vaccine)
initial(IOxGetDie2[,]) <- IOxGetDie2_0[i,j]

gamma_get_ox_die <- user() # rate of progression through requiring oxygen compartment conditional on getting oxygen and dying

deriv(IOxGetDie1[,1]) <- (prob_non_severe_death_treatment[i] * number_get_Ox[i,j]) - gamma_get_ox_die * IOxGetDie1[i,j]
deriv(IOxGetDie1[,2]) <- (prob_non_severe_death_treatment[i] * number_get_Ox[i,j]) - gamma_get_ox_die * IOxGetDie1[i,j] - (gamma_vaccine[j] * IOxGetDie1[i,j])
deriv(IOxGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxGetDie1[i,j-1]) + (prob_non_severe_death_treatment[i] * number_get_Ox[i,j]) - gamma_get_ox_die * IOxGetDie1[i,j] - (gamma_vaccine[j] * IOxGetDie1[i,j])


deriv(IOxGetDie2[,]) <- (gamma_get_ox_die * IOxGetDie1[i,j]) - (gamma_get_ox_die * IOxGetDie2[i,j])
deriv(IOxGetDie2[,2]) <- (gamma_get_ox_die * IOxGetDie1[i,j]) - (gamma_get_ox_die * IOxGetDie2[i,j]) - (gamma_vaccine[j] * IOxGetDie2[i,j])
deriv(IOxGetDie2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxGetDie2[i,j-1]) + (gamma_get_ox_die * IOxGetDie1[i,j]) - (gamma_get_ox_die * IOxGetDie2[i,j]) - (gamma_vaccine[j] * IOxGetDie2[i,j])
################################################################################

### IOxNotGetLive (IOxNotGetLive1 & IOxNotGetLive2): Do not get oxygen, go on to survive #######
dim(IOxNotGetLive1) <- c(N_age, N_vaccine)
dim(IOxNotGetLive2) <- c(N_age, N_vaccine)

IOxNotGetLive1_0[,] <- user()
dim(IOxNotGetLive1_0) <- c(N_age, N_vaccine)
initial(IOxNotGetLive1[,]) <- IOxNotGetLive1_0[i,j]

IOxNotGetLive2_0[,] <- user()
dim(IOxNotGetLive2_0) <- c(N_age, N_vaccine)
initial(IOxNotGetLive2[,]) <- IOxNotGetLive2_0[i,j]

gamma_not_get_ox_survive <- user() # rate of progression through requiring oxygen compartment conditional on not getting oxygen and surviving

deriv(IOxNotGetLive1[,1]) <- ((number_requiring_Ox[i,j] - number_get_Ox[i,j]) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i,j])
deriv(IOxNotGetLive1[,2]) <- ((number_requiring_Ox[i,j] - number_get_Ox[i,j]) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) - (gamma_vaccine[j] * IOxNotGetLive1[i,j])
deriv(IOxNotGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxNotGetLive1[i,j-1]) + ((number_requiring_Ox[i,j] - number_get_Ox[i,j]) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) - (gamma_vaccine[j] * IOxNotGetLive1[i,j])

deriv(IOxNotGetLive2[,1]) <- (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) -  (gamma_not_get_ox_survive * IOxNotGetLive2[i,j])
deriv(IOxNotGetLive2[,2]) <- (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) -  (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) - (gamma_vaccine[j] * IOxNotGetLive2[i,j])
deriv(IOxNotGetLive2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxNotGetLive2[i,j-1]) + (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) -  (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) - (gamma_vaccine[j] * IOxNotGetLive2[i,j])
################################################################################

### IOxNotGetDie (IOxNotGetDie1 & IOxNotGetDie2): Do not get oxygen, go on to die #######
dim(IOxNotGetDie1) <- c(N_age, N_vaccine)
dim(IOxNotGetDie2) <- c(N_age, N_vaccine)

IOxNotGetDie1_0[,] <- user()
dim(IOxNotGetDie1_0) <- c(N_age, N_vaccine)
initial(IOxNotGetDie1[,]) <- IOxNotGetDie1_0[i,j]

IOxNotGetDie2_0[,] <- user()
dim(IOxNotGetDie2_0) <- c(N_age, N_vaccine)
initial(IOxNotGetDie2[,]) <- IOxNotGetDie2_0[i,j]

gamma_not_get_ox_die <- user() # rate of progression through requiring oxygen compartment conditional on not getting oxygen and dying

deriv(IOxNotGetDie1[,1]) <- ((number_requiring_Ox[i,j] - number_get_Ox[i,j]) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i,j])
deriv(IOxNotGetDie1[,2]) <- ((number_requiring_Ox[i,j] - number_get_Ox[i,j]) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_vaccine[j] * IOxNotGetDie1[i,j])
deriv(IOxNotGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxNotGetDie1[i,j-1]) + ((number_requiring_Ox[i,j] - number_get_Ox[i,j]) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_vaccine[j] * IOxNotGetDie1[i,j])

deriv(IOxNotGetDie2[,1]) <- (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_not_get_ox_die * IOxNotGetDie2[i,j])
deriv(IOxNotGetDie2[,2]) <- (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) - (gamma_vaccine[j] * IOxNotGetDie2[i,j])
deriv(IOxNotGetDie2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxNotGetDie2[i,j-1]) + (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) - (gamma_vaccine[j] * IOxNotGetDie2[i,j])
################################################################################

### IMVGetLive (IMVGetLive1 & IMVGetLive2): Get mechanical ventilation, go on to live ########
dim(IMVGetLive1) <- c(N_age, N_vaccine)
dim(IMVGetLive2) <- c(N_age, N_vaccine)

IMVGetLive1_0[,] <- user()
dim(IMVGetLive1_0) <- c(N_age, N_vaccine)
initial(IMVGetLive1[,]) <- IMVGetLive1_0[i,j]

IMVGetLive2_0[,] <- user()
dim(IMVGetLive2_0) <- c(N_age, N_vaccine)
initial(IMVGetLive2[,]) <- IMVGetLive2_0[i,j]

gamma_get_mv_survive <- user() # rate of progression through requiring mechanical ventilation compartment conditional on getting ventilation and surviving

deriv(IMVGetLive1[,1]) <- ((1 - prob_severe_death_treatment[i]) * number_get_IMV[i,j]) - (gamma_get_mv_survive * IMVGetLive1[i,j])
deriv(IMVGetLive1[,2]) <- ((1 - prob_severe_death_treatment[i]) * number_get_IMV[i,j]) - (gamma_get_mv_survive * IMVGetLive1[i,j]) - (gamma_vaccine[j] * IMVGetLive1[i,j])
deriv(IMVGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetLive1[i,j-1]) + ((1 - prob_severe_death_treatment[i]) * number_get_IMV[i,j]) - (gamma_get_mv_survive * IMVGetLive1[i,j]) - (gamma_vaccine[j] * IMVGetLive1[i,j])

deriv(IMVGetLive2[,1]) <- (gamma_get_mv_survive * IMVGetLive1[i,j]) - (gamma_get_mv_survive * IMVGetLive2[i,j])
deriv(IMVGetLive2[,2]) <- (gamma_get_mv_survive * IMVGetLive1[i,j]) - (gamma_get_mv_survive * IMVGetLive2[i,j]) - (gamma_vaccine[j] * IMVGetLive2[i,j])
deriv(IMVGetLive2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetLive2[i,j-1]) + (gamma_get_mv_survive * IMVGetLive1[i,j]) - (gamma_get_mv_survive * IMVGetLive2[i,j]) - (gamma_vaccine[j] * IMVGetLive2[i,j])
################################################################################

### IMVGetDie (IMVGetDie1 & IMVGetDie2): Get mechanical ventilation, go on to die ########
dim(IMVGetDie1) <- c(N_age, N_vaccine)
dim(IMVGetDie2) <- c(N_age, N_vaccine)

IMVGetDie1_0[,] <- user()
dim(IMVGetDie1_0) <- c(N_age, N_vaccine)
initial(IMVGetDie1[,]) <- IMVGetDie1_0[i,j]

IMVGetDie2_0[,] <- user()
dim(IMVGetDie2_0) <- c(N_age, N_vaccine)
initial(IMVGetDie2[,]) <- IMVGetDie2_0[i,j]

gamma_get_mv_die <- user() # rate of progression through requiring mechanical ventilation compartment conditional on getting ventilation and dying

deriv(IMVGetDie1[,1]) <- (prob_severe_death_treatment[i] * number_get_IMV[i,j]) - (gamma_get_mv_die * IMVGetDie1[i,j])
deriv(IMVGetDie1[,2]) <- (prob_severe_death_treatment[i] * number_get_IMV[i,j]) - (gamma_get_mv_die * IMVGetDie1[i,j]) - (gamma_vaccine[j] * IMVGetDie1[i,j])
deriv(IMVGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetDie1[i,j-1]) + (prob_severe_death_treatment[i] * number_get_IMV[i,j]) - (gamma_get_mv_die * IMVGetDie1[i,j]) - (gamma_vaccine[j] * IMVGetDie1[i,j])

deriv(IMVGetDie2[,1]) <- (gamma_get_mv_die * IMVGetDie1[i,j]) - (gamma_get_mv_die * IMVGetDie2[i,j])
deriv(IMVGetDie2[,2]) <- (gamma_get_mv_die * IMVGetDie1[i,j]) - (gamma_get_mv_die * IMVGetDie2[i,j]) - (gamma_vaccine[j] * IMVGetDie2[i,j])
deriv(IMVGetDie2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetDie2[i,j-1]) + (gamma_get_mv_die * IMVGetDie1[i,j]) - (gamma_get_mv_die * IMVGetDie2[i,j]) - (gamma_vaccine[j] * IMVGetDie2[i,j])
################################################################################

### IMVNotGetLive (IMVNotGetLive1 & IMVNotGetLive2): Do no get mechanical ventilation, go on to live ########
dim(IMVNotGetLive1) <- c(N_age, N_vaccine)
dim(IMVNotGetLive2) <- c(N_age, N_vaccine)

IMVNotGetLive1_0[,] <- user()
dim(IMVNotGetLive1_0) <- c(N_age, N_vaccine)
initial(IMVNotGetLive1[,]) <- IMVNotGetLive1_0[i,j]

IMVNotGetLive2_0[,] <- user()
dim(IMVNotGetLive2_0) <- c(N_age, N_vaccine)
initial(IMVNotGetLive2[,]) <- IMVNotGetLive2_0[i,j]

gamma_not_get_mv_survive <- user() # rate of progression through requiring mechanical ventilation compartment conditional on not getting ventilation and surviving

deriv(IMVNotGetLive1[,1]) <- ((number_requiring_IMV[i,j] - number_get_IMV[i,j]) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i,j])
deriv(IMVNotGetLive1[,2]) <- ((number_requiring_IMV[i,j] - number_get_IMV[i,j]) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_vaccine[j] * IMVNotGetLive1[i,j])
deriv(IMVNotGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVNotGetLive1[i,j-1]) + ((number_requiring_IMV[i,j] - number_get_IMV[i,j]) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_vaccine[j] * IMVNotGetLive1[i,j])

deriv(IMVNotGetLive2[,1]) <- (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i,j])
deriv(IMVNotGetLive2[,2]) <- (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_vaccine[j] * IMVNotGetLive2[i,j])
deriv(IMVNotGetLive2[,3:N_vaccine]) <-  - (gamma_vaccine[j-1] * IMVNotGetLive2[i,j-1]) + (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_vaccine[j] * IMVNotGetLive2[i,j])
################################################################################

### IMVNotGetDie (IMVNotGetDie1 & IMVNotGetDie2): Do no get mechanical ventilation, go on to die ########
dim(IMVNotGetDie1) <- c(N_age, N_vaccine)
dim(IMVNotGetDie2) <- c(N_age, N_vaccine)

IMVNotGetDie1_0[,] <- user()
dim(IMVNotGetDie1_0) <- c(N_age, N_vaccine)
initial(IMVNotGetDie1[,]) <- IMVNotGetDie1_0[i,j]

IMVNotGetDie2_0[,] <- user()
dim(IMVNotGetDie2_0) <- c(N_age, N_vaccine)
initial(IMVNotGetDie2[,]) <- IMVNotGetDie2_0[i,j]

gamma_not_get_mv_die <- user() # rate of progression through requiring mechanical ventilation compartment conditional on not getting ventilation and dying

deriv(IMVNotGetDie1[,1]) <- ((number_requiring_IMV[i,j] - number_get_IMV[i,j]) * prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i,j])
deriv(IMVNotGetDie1[,2]) <- ((number_requiring_IMV[i,j] - number_get_IMV[i,j]) * prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_vaccine[j] * IMVNotGetDie1[i,j])
deriv(IMVNotGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVNotGetDie1[i,j-1]) + ((number_requiring_IMV[i,j] - number_get_IMV[i,j]) * prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_vaccine[j] * IMVNotGetDie1[i,j])

deriv(IMVNotGetDie2[,1]) <- (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_not_get_mv_die * IMVNotGetDie2[i,j])
deriv(IMVNotGetDie2[,2]) <- (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_not_get_mv_die * IMVNotGetDie2[i,j]) - (gamma_vaccine[j] * IMVNotGetDie2[i,j])
deriv(IMVNotGetDie2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVNotGetDie2[i,j-1]) + (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_not_get_mv_die * IMVNotGetDie2[i,j]) - (gamma_vaccine[j] * IMVNotGetDie2[i,j])
################################################################################

### IRec (IRec1 & IRec2): Recovering from ICU ##################################
dim(IRec1) <- c(N_age, N_vaccine)
dim(IRec2) <- c(N_age, N_vaccine)
dim(IRec) <- N_age

IRec1_0[,] <- user()
dim(IRec1_0) <- c(N_age, N_vaccine)
initial(IRec1[,]) <- IRec1_0[i,j]

IRec2_0[,] <- user()
dim(IRec2_0) <- c(N_age, N_vaccine)
initial(IRec2[,]) <- IRec2_0[i,j]

gamma_rec <- user() # rate of progression through post-ICU recovery compartment

deriv(IRec1[,1]) <- (gamma_get_mv_survive * IMVGetLive2[i,j]) - (gamma_rec * IRec1[i,j])
deriv(IRec1[,2]) <- (gamma_get_mv_survive * IMVGetLive2[i,j]) - (gamma_rec * IRec1[i,j]) - (gamma_vaccine[j] * IRec1[i,j])
deriv(IRec1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IRec1[i,j-1]) + (gamma_get_mv_survive * IMVGetLive2[i,j]) - (gamma_rec * IRec1[i,j]) - (gamma_vaccine[j] * IRec1[i,j])

deriv(IRec2[,1]) <- (gamma_rec * IRec1[i,j]) - (gamma_rec * IRec2[i,j])
deriv(IRec2[,2]) <- (gamma_rec * IRec1[i,j]) - (gamma_rec * IRec2[i,j]) - (gamma_vaccine[j] * IRec2[i,j])
deriv(IRec2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IRec2[i,j-1]) + (gamma_rec * IRec1[i,j]) - (gamma_rec * IRec2[i,j]) - (gamma_vaccine[j] * IRec2[i,j])

output(IRec[]) <- sum(IRec1[i,]) + sum(IRec2[i,])
################################################################################

### D: Dead ####################################################################
dim(D) <- c(N_age, N_vaccine)

D_0[,] <- user()
dim(D_0) <- c(N_age, N_vaccine)
initial(D[,]) <- D_0[i,j]

deriv(D[,1]) <- (gamma_get_ox_die * IOxGetDie2[i,j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) + (gamma_get_mv_die * IMVGetDie2[i,j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i,j])
deriv(D[,2]) <- (gamma_get_ox_die * IOxGetDie2[i,j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) + (gamma_get_mv_die * IMVGetDie2[i,j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i,j]) - (gamma_vaccine[j] * D[i,j])
deriv(D[,3:N_vaccine]) <- (gamma_vaccine[j-1] * D[i,j-1]) + (gamma_get_ox_die * IOxGetDie2[i,j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) + (gamma_get_mv_die * IMVGetDie2[i,j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i,j]) - (gamma_vaccine[j] * D[i,j])
################################################################################

################################################################################
### Vaccination capacity #######################################################
################################################################################
# Vaccination
vaccination_target[] <- user() # 0/1 index of targeted age groups
dim(vaccination_target) <- N_age

vaccine_efficacy_infection[,] <- user() # Reduction in lambda from vaccination by age and vaccination status
dim(vaccine_efficacy_infection) <- c(N_age, N_vaccine)

gamma_vaccine[] <- user() # Vector of vaccine progression parameters by vaccination status (First = 0 as handled separately as time-varying vaccination rate, Last = 0 as no progression from "previously vaccinated group)
dim(gamma_vaccine) <- N_vaccine

# Interpolation of vaccination rate over time
mv <- interpolate(tt_vaccine, max_vaccine, "constant")
tt_vaccine[] <- user()
max_vaccine[] <- user()
dim(tt_vaccine) <- user()
dim(max_vaccine) <- length(tt_vaccine)
vr_temp[] <- S[i,1] * vaccination_target[i] + E1[i,1] * vaccination_target[i] + E2[i,1] * vaccination_target[i] + R1[i,1] * vaccination_target[i] + R2[i,1] * vaccination_target[i]
dim(vr_temp) <- N_age
# Catch so max vaccination rate does not -> infinity as vaccine-eligible population -> 0
vr_den <- if(sum(vr_temp) <= mv) mv else sum(vr_temp)
vr <- mv / vr_den # Vaccination rate to achieve capacity given number in vaccine-eligible population
################################################################################
################################################################################

################################################################################
### Hospital and ICU capacity ##################################################
################################################################################
## Interpolation for Hospital and ICU Capacity
hosp_bed_capacity <- interpolate(tt_hosp_beds, hosp_beds, "constant")
tt_hosp_beds[] <- user()
hosp_beds[] <- user()
dim(tt_hosp_beds) <- user()
dim(hosp_beds) <- length(tt_hosp_beds)

ICU_bed_capacity <- interpolate(tt_ICU_beds, ICU_beds, "constant")
tt_ICU_beds[] <- user()
ICU_beds[] <- user()
dim(tt_ICU_beds) <- user()
dim(ICU_beds) <- length(tt_ICU_beds)

prob_hosp[,] <- user() # probability of requiring hospitalisation by age and vaccination status
dim(prob_hosp) <- c(N_age, N_vaccine)

prob_severe[] <- user() # probability of severe disease (requiring mechanical ventilation) by age
dim(prob_severe) <- N_age

prob_non_severe_death_treatment[] <- user() # probability of dying from non-severe disease (i.e. requiring oxygen but not mechanical ventilation) by age given you receive appropriate treatment (proxy here is whether a general hospital bed is available)
dim(prob_non_severe_death_treatment) <- N_age

prob_non_severe_death_no_treatment[] <- user() # probability of dying from non-severe disease (i.e. requiring oxygen but not mechanical ventilation) by age given you do NOT receive appropriate treatment (proxy here is whether a general hospital bed is available)
dim(prob_non_severe_death_no_treatment) <- N_age

prob_severe_death_treatment[] <- user() # probability of dying from severe disease (i.e. requiring mechanical ventilation) by age given you receive appropriate treatment (proxy here is whether an ICU bed is available)
dim(prob_severe_death_treatment) <- N_age

prob_severe_death_no_treatment[] <- user() # probability of dying from severe disease (i.e. requiring mechanical ventilation) by age given you do NOT receive appropriate treatment (proxy here is whether an ICU bed is available)
dim(prob_severe_death_no_treatment) <- N_age

p_dist[,] <- user() # distributing infections in given age/vaccine class to available hosp/ICU beds (make all equal to make it random and not related to age)
dim(p_dist) <- c(N_age, N_vaccine)

# Infections Requiring Oxygen (a general Hosptial Bed)
hosp_occ <- sum(IOxGetLive1) + sum(IOxGetLive2) + sum(IOxGetDie1) + sum(IOxGetDie2) + sum(IRec1) + sum(IRec2) # Summing number of infections in compartments that use general hospital beds
current_free_hosp <- hosp_bed_capacity - hosp_occ + gamma_get_ox_die*sum(IOxGetDie2) + gamma_get_ox_survive * sum(IOxGetLive2) + gamma_rec * sum(IRec2) - gamma_get_mv_survive * sum(IMVGetLive2)

number_requiring_Ox[,] <- gamma_ICase * ICase2[i,j] * (1 - prob_severe[i]) # NOTE THIS IS DIFF IN SYNTAX FROM STOCHSTIC VERSION WHERE WE SUBTRACT THE NUMBER GETTING IMV - MIGHT BE BETTER FROM A ROUNDING ERROR PERSPECITVE
dim(number_requiring_Ox) <- c(N_age, N_vaccine)

total_number_requiring_ox <- sum(number_requiring_Ox)
total_number_get_hosp <- if (current_free_hosp <= 0) 0 else (if(current_free_hosp - total_number_requiring_ox >= 0) total_number_requiring_ox else(current_free_hosp)) # Working out the number of new hospital bed requiring infections that get a bed

Ox_dist_weighting[,] <- number_requiring_Ox[i,j] * p_dist[i,j]
dim(Ox_dist_weighting) <- c(N_age, N_vaccine)

number_get_Ox[,] <- if (total_number_requiring_ox == 0) 0 else Ox_dist_weighting[i,j]/sum(Ox_dist_weighting) * total_number_get_hosp
dim(number_get_Ox) <- c(N_age, N_vaccine)

# Infections Requiring Mechanical Ventilation (an ICU Bed)
ICU_occ <- sum(IMVGetLive1) + sum(IMVGetLive2) + sum(IMVGetDie1) + sum(IMVGetDie2) # Summing number of infections in compartments that use ICU beds
current_free_ICUs <- ICU_bed_capacity - ICU_occ + gamma_get_mv_survive *sum(IMVGetLive2) + gamma_get_mv_die *sum(IMVGetDie2)

number_requiring_IMV[,] <- gamma_ICase * ICase2[i,j] * prob_severe[i]
dim(number_requiring_IMV) <- c(N_age, N_vaccine)

total_number_requiring_IMV <- sum(number_requiring_IMV)
total_number_get_IMV <- if(current_free_ICUs <= 0) 0 else(if(current_free_ICUs - total_number_requiring_IMV >= 0) total_number_requiring_IMV else(current_free_ICUs)) # Working out the number of new ICU requiring infections that get a bed

IMV_dist_weighting[,] <- number_requiring_IMV[i,j] * p_dist[i,j]
dim(IMV_dist_weighting) <- c(N_age, N_vaccine)

number_get_IMV[,] <- if (total_number_requiring_IMV == 0) 0 else IMV_dist_weighting[i,j]/sum(IMV_dist_weighting) * total_number_get_IMV
dim(number_get_IMV) <- c(N_age, N_vaccine)
################################################################################
################################################################################

################################################################################
### FOI and contact matrix #####################################################
################################################################################
# Generating Force of Infection
m[, ] <- interpolate(tt_matrix, mix_mat_set, "constant")
dim(m) <- c(N_age, N_age)
tt_matrix[] <- user()
mix_mat_set[, ,] <- user()
dim(tt_matrix) <- user()
dim(mix_mat_set) <- c(length(tt_matrix), N_age, N_age)

# Interpolation for beta
beta <- interpolate(tt_beta, beta_set, "constant")
tt_beta[] <- user()
beta_set[] <- user()
dim(tt_beta) <- user()
dim(beta_set) <- length(tt_beta)

# Generating Force of Infection
temp[] <- sum(IMild[i,]) + sum(ICase1[i,]) + sum(ICase2[i,])
dim(temp) <- N_age

s_ij[,] <- m[i, j] * temp[j]
dim(s_ij) <- c(N_age, N_age)

lambda[] <- beta * sum(s_ij[i, ])
dim(lambda) <- N_age
################################################################################
################################################################################

################################################################################
### Output #####################################################################
################################################################################
# Hospital occupancy and demand
output(hospital_occupancy[]) <- sum(IOxGetLive1[i,]) + sum(IOxGetLive2[i,]) + sum(IOxGetDie1[i,]) + sum(IOxGetDie2[i,]) + sum(IRec1[i,]) + sum(IRec2[i,])
dim(hospital_occupancy) <- N_age

output(ICU_occupancy[]) <- sum(IMVGetLive1[i,]) + sum(IMVGetLive2[i,]) + sum(IMVGetDie1[i,]) + sum(IMVGetDie2[i,])
dim(ICU_occupancy) <- N_age

output(hospital_demand[]) <- sum(IOxGetLive1[i,]) + sum(IOxGetLive2[i,]) + sum(IOxGetDie1[i,]) + sum(IOxGetDie2[i,]) + sum(IRec1[i,]) + sum(IRec2[i,]) + sum(IOxNotGetLive1[i,]) + sum(IOxNotGetLive2[i,]) + sum(IOxNotGetDie1[i,]) + sum(IOxNotGetDie2[i,])
dim(hospital_demand) <- N_age

output(ICU_demand[]) <- sum(IMVGetLive1[i,]) +  sum(IMVGetLive2[i,]) + sum(IMVGetDie1[i,]) + sum(IMVGetDie2[i,]) + sum(IMVNotGetLive1[i,]) + sum(IMVNotGetLive2[i,]) + sum(IMVNotGetDie1[i,]) + sum(IMVNotGetDie2[i,])
dim(ICU_demand) <- N_age

# Number in hospital or ICU compartments
output(IICU[]) <- sum(IMVGetLive1[i,]) + sum(IMVGetLive2[i,]) + sum(IMVGetDie1[i,]) + sum(IMVGetDie2[i,]) + sum(IMVNotGetLive1[i,]) + sum(IMVNotGetLive2[i,]) + sum(IMVNotGetDie1[i,]) + sum(IMVNotGetDie2[i,])
dim(IICU) <- N_age

output(IHospital[]) <- sum(IOxGetLive1[i,]) + sum(IOxGetLive2[i,]) + sum(IOxGetDie1[i,]) + sum(IOxGetDie2[i,]) + sum(IOxNotGetLive1[i,]) + sum(IOxNotGetLive2[i,]) + sum(IOxNotGetDie1[i,]) + sum(IOxNotGetDie2[i,])
dim(IHospital) <- N_age

# Deaths
Dlag[,] <- delay(D[i,j], dt)
dim(Dlag) <- c(N_age, N_vaccine)
output(deaths[]) <- sum(D[i,]) - sum(Dlag[i,])
dim(deaths) <- N_age

# Infections
deriv(I[,]) <- (lambda[i] * vaccine_efficacy_infection[i,j] * S[i,j])
dim(I) <- c(N_age, N_vaccine)
initial(I[,]) <- 0

Ilag[,] <- delay(I[i,j], dt)
dim(Ilag) <- c(N_age, N_vaccine)
output(infections[]) <- sum(I[i,]) - sum(Ilag[i,])
dim(infections) <- N_age

# Vaccinations
deriv(V[]) <- (vr * vaccination_target[i] * S[i,1]) + (vr * vaccination_target[i] * E1[i,1]) + (vr * vaccination_target[i] * E2[i,1]) +
  (vr * vaccination_target[i] * R2[i,1]) + (vr * vaccination_target[i] * R2[i,1])
dim(V) <- N_age
initial(V[]) <- 0

Vlag[] <- delay(V[i], dt)
dim(Vlag) <- N_age
output(vaccines[]) <- V[i] - Vlag[i]
dim(vaccines) <- N_age

# Population size
output(N[]) <- sum(S[i,]) + sum(E1[i,]) + sum(E2[i,]) + sum(IMild[i,]) + sum(ICase1[i,]) + sum(ICase2[i,]) +
  sum(IMVGetLive1[i,]) + sum(IMVGetLive2[i,]) +
  sum(IMVGetDie1[i,]) + sum(IMVGetDie2[i,]) + sum(IMVNotGetLive1[i,]) + sum(IMVNotGetLive2[i,]) + sum(IMVNotGetDie1[i,]) + sum(IMVNotGetDie2[i,]) +
  sum(IOxGetLive1[i,]) + sum(IOxGetLive2[i,]) + sum(IOxGetDie1[i,]) + sum(IOxGetDie2[i,]) + sum(IOxNotGetLive1[i,]) + sum(IOxNotGetLive2[i,]) +
  sum(IOxNotGetDie1[i,]) + sum(IOxNotGetDie2[i,]) +
  sum(IRec1[i,]) + sum(IRec2[i,]) +
  sum(R1[i,]) + sum(R2[i,]) + sum(D[i,])
dim(N) <- N_age
################################################################################