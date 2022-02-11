################################################################################
### Vaccine model: deterministic ###############################################
################################################################################

### Initial setup ##############################################################
N_age <- user() # Number of age groups
N_vaccine <- user() # Number of vaccine groups
################################################################################

## Time output to line up with squire fitting infrastructure
time <- t
output(time) <- TRUE


### S: susceptible #############################################################
dim(S) <- c(N_age, N_vaccine)

S_0[,] <- user()
dim(S_0) <- c(N_age, N_vaccine)
initial(S[,]) <- S_0[i,j]

deriv(S[,1]) <- (gamma_R_t * R2[i,j]) - (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (vr * vaccination_target[i] * S[i,j])
deriv(S[,2]) <- (vr * vaccination_target[i] * S[i,j-1]) + (gamma_R_t * R2[i,j]) - (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_vaccine[j] * S[i,j])
deriv(S[,3:N_vaccine]) <- (gamma_vaccine[j-1] * S[i,j-1]) + (gamma_R_t * R2[i,j]) - (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_vaccine[j] * S[i,j])
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

deriv(E1[,1]) <- (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_E * E1[i,j]) - (vr * vaccination_target[i] * E1[i,j])
deriv(E1[,2]) <- (vr * vaccination_target[i] * E1[i,j-1]) + (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_E * E1[i,j]) - (gamma_vaccine[j] * E1[i,j])
deriv(E1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * E1[i,j-1]) + (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j]) - (gamma_E * E1[i,j]) - (gamma_vaccine[j] * E1[i,j])

deriv(E2[,1]) <- (gamma_E * E1[i,j]) - (gamma_E * E2[i,j]) - (vr * vaccination_target[i] * E2[i,j])
deriv(E2[,2]) <- (vr * vaccination_target[i] * E2[i,j-1]) + (gamma_E * E1[i,j]) - (gamma_E * E2[i,j]) - (gamma_vaccine[j] * E2[i,j])
deriv(E2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * E2[i,j-1]) + (gamma_E * E1[i,j]) - (gamma_E * E2[i,j]) - (gamma_vaccine[j] * E2[i,j])

output(E[]) <- sum(E1[i,]) + sum(E2[i,])
dim(E) <- N_age
################################################################################

### IMild: Unhospitalised infection ############################################
dim(IMild) <- c(N_age, N_vaccine)

IMild_0[,] <- user()
dim(IMild_0) <- c(N_age, N_vaccine)
initial(IMild[,]) <- IMild_0[i,j]

gamma_IMild <- user() # rate of progression from mild infection to recovery

deriv(IMild[,1]) <- (gamma_E * E2[i,j] * (1 - prob_hosp_t_mult[i,j])) - (gamma_IMild * IMild[i,j])
deriv(IMild[,2]) <- (gamma_E * E2[i,j] * (1 - prob_hosp_t_mult[i,j])) - (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * IMild[i,j])
deriv(IMild[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMild[i,j-1]) + (gamma_E * E2[i,j] * (1 - prob_hosp_t_mult[i,j])) - (gamma_IMild * IMild[i,j]) - (gamma_vaccine[j] * IMild[i,j])

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

# Interpolation for dur_R
gamma_R_t <- interpolate(tt_dur_R, gamma_R, "constant")
tt_dur_R[] <- user()
dim(tt_dur_R) <- user()
dim(gamma_R) <- length(tt_dur_R)
gamma_R[] <- user() # rate of progression through recovered compartment (loss of naturally acquired immunity)

deriv(R1[,1]) <- (gamma_rec * IRec2[i,j]) + (gamma_IMild * IMild[i,j]) + (gamma_get_ox_survive_t * IOxGetLive2[i,j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_R_t * R1[i,j]) - (vr * vaccination_target[i] * R1[i,j])
deriv(R1[,2]) <- (vr * vaccination_target[i] * R1[i,j-1]) + (gamma_rec * IRec2[i,j]) + (gamma_IMild * IMild[i,j]) + (gamma_get_ox_survive_t * IOxGetLive2[i,j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_R_t * R1[i,j]) - (gamma_vaccine[j] * R1[i,j])
deriv(R1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * R1[i,j-1]) + (gamma_rec * IRec2[i,j]) + (gamma_IMild * IMild[i,j]) + (gamma_get_ox_survive_t * IOxGetLive2[i,j]) + (gamma_not_get_ox_survive * IOxNotGetLive2[i,j]) + (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_R_t * R1[i,j]) - (gamma_vaccine[j] * R1[i,j])

deriv(R2[,1]) <- (gamma_R_t * R1[i,j]) - (gamma_R_t * R2[i,j]) - (vr * vaccination_target[i] * R2[i,j])
deriv(R2[,2]) <- (vr * vaccination_target[i] * R2[i,j-1]) + (gamma_R_t * R1[i,j]) - (gamma_R_t * R2[i,j]) - (gamma_vaccine[j] * R2[i,j])
deriv(R2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * R2[i,j-1]) + (gamma_R_t * R1[i,j]) - (gamma_R_t * R2[i,j]) - (gamma_vaccine[j] * R2[i,j])

output(R[]) <- sum(R1[i,]) + sum(R2[i,])
dim(R) <- N_age
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

deriv(ICase1[,1]) <- (gamma_E * E2[i,j] * prob_hosp_t_mult[i,j]) - (gamma_ICase * ICase1[i,j])
deriv(ICase1[,2]) <- (gamma_E * E2[i,j] * prob_hosp_t_mult[i,j]) - (gamma_ICase * ICase1[i,j]) - (gamma_vaccine[j] * ICase1[i,j])
deriv(ICase1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * ICase1[i,j-1]) + (gamma_E * E2[i,j] * prob_hosp_t_mult[i,j]) - (gamma_ICase * ICase1[i,j]) - (gamma_vaccine[j] * ICase1[i,j])

deriv(ICase2[,1]) <- (gamma_ICase * ICase1[i,j]) - (gamma_ICase * ICase2[i,j])
deriv(ICase2[,2]) <- (gamma_ICase * ICase1[i,j]) - (gamma_ICase * ICase2[i,j]) - (gamma_vaccine[j] * ICase2[i,j])
deriv(ICase2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * ICase2[i,j-1]) + (gamma_ICase * ICase1[i,j]) - (gamma_ICase * ICase2[i,j]) - (gamma_vaccine[j] * ICase2[i,j])

output(ICase[]) <- sum(ICase1[i,]) + sum(ICase2[i,])
dim(ICase) <- N_age
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

gamma_get_ox_survive[] <- user() # rate of progression through requiring oxygen compartment conditional on getting oxygen and surviving
dim(gamma_get_ox_survive) <- user()
tt_dur_get_ox_survive[] <- user()
dim(tt_dur_get_ox_survive) <- length(gamma_get_ox_survive)
gamma_get_ox_survive_t <- interpolate(tt_dur_get_ox_survive, gamma_get_ox_survive, "constant")

deriv(IOxGetLive1[,1]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * p_oxygen * (1 - prob_non_severe_death_treatment[i])) - (gamma_get_ox_survive_t * IOxGetLive1[i,j])
deriv(IOxGetLive1[,2]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * p_oxygen * (1 - prob_non_severe_death_treatment[i])) - (gamma_get_ox_survive_t * IOxGetLive1[i,j]) - (gamma_vaccine[j] * IOxGetLive1[i,j])
deriv(IOxGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxGetLive1[i,j-1]) + (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * p_oxygen * (1 - prob_non_severe_death_treatment[i])) - (gamma_get_ox_survive_t * IOxGetLive1[i,j]) - (gamma_vaccine[j] * IOxGetLive1[i,j])

deriv(IOxGetLive2[,1]) <- (gamma_get_ox_survive_t * IOxGetLive1[i,j]) - (gamma_get_ox_survive_t * IOxGetLive2[i,j])
deriv(IOxGetLive2[,2]) <- (gamma_get_ox_survive_t * IOxGetLive1[i,j]) - (gamma_get_ox_survive_t * IOxGetLive2[i,j]) - (gamma_vaccine[j] *  IOxGetLive2[i,j])
deriv(IOxGetLive2[,3:N_vaccine]) <- (gamma_vaccine[j-1] *  IOxGetLive2[i,j-1]) + (gamma_get_ox_survive_t * IOxGetLive1[i,j]) -  (gamma_get_ox_survive_t * IOxGetLive2[i,j]) - (gamma_vaccine[j] *  IOxGetLive2[i,j])
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

gamma_get_ox_die[] <- user() # rate of progression through requiring oxygen compartment conditional on getting oxygen and dying
dim(gamma_get_ox_die) <- user()
tt_dur_get_ox_die[] <- user()
dim(tt_dur_get_ox_die) <- length(gamma_get_ox_die)
gamma_get_ox_die_t <- interpolate(tt_dur_get_ox_die, gamma_get_ox_die, "constant")

deriv(IOxGetDie1[,1]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * p_oxygen * prob_non_severe_death_treatment[i]) - gamma_get_ox_die_t * IOxGetDie1[i,j]
deriv(IOxGetDie1[,2]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * p_oxygen * prob_non_severe_death_treatment[i])- gamma_get_ox_die_t * IOxGetDie1[i,j] - (gamma_vaccine[j] * IOxGetDie1[i,j])
deriv(IOxGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxGetDie1[i,j-1]) + (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * p_oxygen * prob_non_severe_death_treatment[i]) - gamma_get_ox_die_t * IOxGetDie1[i,j] - (gamma_vaccine[j] * IOxGetDie1[i,j])

deriv(IOxGetDie2[,]) <- (gamma_get_ox_die_t * IOxGetDie1[i,j]) - (gamma_get_ox_die_t * IOxGetDie2[i,j])
deriv(IOxGetDie2[,2]) <- (gamma_get_ox_die_t * IOxGetDie1[i,j]) - (gamma_get_ox_die_t * IOxGetDie2[i,j]) - (gamma_vaccine[j] * IOxGetDie2[i,j])
deriv(IOxGetDie2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxGetDie2[i,j-1]) + (gamma_get_ox_die_t * IOxGetDie1[i,j]) - (gamma_get_ox_die_t * IOxGetDie2[i,j]) - (gamma_vaccine[j] * IOxGetDie2[i,j])
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

deriv(IOxNotGetLive1[,1]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * (1-p_oxygen) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i,j])
deriv(IOxNotGetLive1[,2]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * (1-p_oxygen) * (1 - prob_non_severe_death_no_treatment[i]))- (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) - (gamma_vaccine[j] * IOxNotGetLive1[i,j])
deriv(IOxNotGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxNotGetLive1[i,j-1]) + (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * (1-p_oxygen) * (1 - prob_non_severe_death_no_treatment[i])) - (gamma_not_get_ox_survive * IOxNotGetLive1[i,j]) - (gamma_vaccine[j] * IOxNotGetLive1[i,j])

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

deriv(IOxNotGetDie1[,1]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * (1-p_oxygen) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i,j])
deriv(IOxNotGetDie1[,2]) <- (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * (1-p_oxygen) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_vaccine[j] * IOxNotGetDie1[i,j])
deriv(IOxNotGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IOxNotGetDie1[i,j-1]) + (gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i]) * (1-p_oxygen) * prob_non_severe_death_no_treatment[i]) - (gamma_not_get_ox_die * IOxNotGetDie1[i,j]) - (gamma_vaccine[j] * IOxNotGetDie1[i,j])

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

gamma_get_mv_survive[] <- user() # rate of progression through requiring mechanical ventilation compartment conditional on getting ventilation and surviving
dim(gamma_get_mv_survive) <- user()
tt_dur_get_mv_survive[] <- user()
dim(tt_dur_get_mv_survive) <- length(gamma_get_mv_survive)
gamma_get_mv_survive_t <- interpolate(tt_dur_get_mv_survive, gamma_get_mv_survive, "constant")

deriv(IMVGetLive1[,1]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * p_ventilation * (1 - prob_severe_death_treatment[i])) - (gamma_get_mv_survive_t * IMVGetLive1[i,j])
deriv(IMVGetLive1[,2]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * p_ventilation * (1 - prob_severe_death_treatment[i])) - (gamma_get_mv_survive_t * IMVGetLive1[i,j]) - (gamma_vaccine[j] * IMVGetLive1[i,j])
deriv(IMVGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetLive1[i,j-1]) + (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * p_ventilation * (1 - prob_severe_death_treatment[i])) - (gamma_get_mv_survive_t * IMVGetLive1[i,j]) - (gamma_vaccine[j] * IMVGetLive1[i,j])

deriv(IMVGetLive2[,1]) <- (gamma_get_mv_survive_t * IMVGetLive1[i,j]) - (gamma_get_mv_survive_t * IMVGetLive2[i,j])
deriv(IMVGetLive2[,2]) <- (gamma_get_mv_survive_t * IMVGetLive1[i,j]) - (gamma_get_mv_survive_t * IMVGetLive2[i,j]) - (gamma_vaccine[j] * IMVGetLive2[i,j])
deriv(IMVGetLive2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetLive2[i,j-1]) + (gamma_get_mv_survive_t * IMVGetLive1[i,j]) - (gamma_get_mv_survive_t * IMVGetLive2[i,j]) - (gamma_vaccine[j] * IMVGetLive2[i,j])
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

gamma_get_mv_die[] <- user() # rate of progression through requiring mechanical ventilation compartment conditional on getting ventilation and dying
dim(gamma_get_mv_die) <- user()
tt_dur_get_mv_die[] <- user()
dim(tt_dur_get_mv_die) <- length(gamma_get_mv_die)
gamma_get_mv_die_t <- interpolate(tt_dur_get_mv_die, gamma_get_mv_die, "constant")

deriv(IMVGetDie1[,1]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * p_ventilation * prob_severe_death_treatment[i]) - (gamma_get_mv_die_t * IMVGetDie1[i,j])
deriv(IMVGetDie1[,2]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * p_ventilation * prob_severe_death_treatment[i]) - (gamma_get_mv_die_t * IMVGetDie1[i,j]) - (gamma_vaccine[j] * IMVGetDie1[i,j])
deriv(IMVGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetDie1[i,j-1]) + (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * p_ventilation * prob_severe_death_treatment[i]) - (gamma_get_mv_die_t * IMVGetDie1[i,j]) - (gamma_vaccine[j] * IMVGetDie1[i,j])

deriv(IMVGetDie2[,1]) <- (gamma_get_mv_die_t * IMVGetDie1[i,j]) - (gamma_get_mv_die_t * IMVGetDie2[i,j])
deriv(IMVGetDie2[,2]) <- (gamma_get_mv_die_t * IMVGetDie1[i,j]) - (gamma_get_mv_die_t * IMVGetDie2[i,j]) - (gamma_vaccine[j] * IMVGetDie2[i,j])
deriv(IMVGetDie2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVGetDie2[i,j-1]) + (gamma_get_mv_die_t * IMVGetDie1[i,j]) - (gamma_get_mv_die_t * IMVGetDie2[i,j]) - (gamma_vaccine[j] * IMVGetDie2[i,j])
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

deriv(IMVNotGetLive1[,1]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * (1 - p_ventilation) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i,j])
deriv(IMVNotGetLive1[,2]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * (1 - p_ventilation) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_vaccine[j] * IMVNotGetLive1[i,j])
deriv(IMVNotGetLive1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVNotGetLive1[i,j-1]) + (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * (1 - p_ventilation) * (1 - prob_severe_death_no_treatment[i])) - (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_vaccine[j] * IMVNotGetLive1[i,j])

deriv(IMVNotGetLive2[,1]) <- (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i,j])
deriv(IMVNotGetLive2[,2]) <- (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_vaccine[j] * IMVNotGetLive2[i,j])
deriv(IMVNotGetLive2[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVNotGetLive2[i,j-1]) + (gamma_not_get_mv_survive * IMVNotGetLive1[i,j]) - (gamma_not_get_mv_survive * IMVNotGetLive2[i,j]) - (gamma_vaccine[j] * IMVNotGetLive2[i,j])
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

deriv(IMVNotGetDie1[,1]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * (1 - p_ventilation) *  prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i,j])
deriv(IMVNotGetDie1[,2]) <- (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * (1 - p_ventilation) *  prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_vaccine[j] * IMVNotGetDie1[i,j])
deriv(IMVNotGetDie1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IMVNotGetDie1[i,j-1]) + (gamma_ICase * ICase2[i,j] * prob_severe_multi[i] * (1 - p_ventilation) *  prob_severe_death_no_treatment[i]) - (gamma_not_get_mv_die * IMVNotGetDie1[i,j]) - (gamma_vaccine[j] * IMVNotGetDie1[i,j])

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

deriv(IRec1[,1]) <- (gamma_get_mv_survive_t * IMVGetLive2[i,j]) - (gamma_rec * IRec1[i,j])
deriv(IRec1[,2]) <- (gamma_get_mv_survive_t * IMVGetLive2[i,j]) - (gamma_rec * IRec1[i,j]) - (gamma_vaccine[j] * IRec1[i,j])
deriv(IRec1[,3:N_vaccine]) <- (gamma_vaccine[j-1] * IRec1[i,j-1]) + (gamma_get_mv_survive_t * IMVGetLive2[i,j]) - (gamma_rec * IRec1[i,j]) - (gamma_vaccine[j] * IRec1[i,j])

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

deriv(D[,1]) <- (gamma_get_ox_die_t * IOxGetDie2[i,j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) + (gamma_get_mv_die_t * IMVGetDie2[i,j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i,j])
deriv(D[,2]) <- (gamma_get_ox_die_t * IOxGetDie2[i,j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) + (gamma_get_mv_die_t * IMVGetDie2[i,j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i,j]) - (gamma_vaccine[j] * D[i,j])
deriv(D[,3:N_vaccine]) <- (gamma_vaccine[j-1] * D[i,j-1]) + (gamma_get_ox_die_t * IOxGetDie2[i,j]) + (gamma_not_get_ox_die * IOxNotGetDie2[i,j]) + (gamma_get_mv_die_t * IMVGetDie2[i,j]) + (gamma_not_get_mv_die * IMVNotGetDie2[i,j]) - (gamma_vaccine[j] * D[i,j])
################################################################################

################################################################################
### Vaccination capacity #######################################################
################################################################################
# Vaccination
# Vaccine prioritisation coverage matrix
N_prioritisation_steps <- user()
vaccine_coverage_mat[,] <- user()
dim(vaccine_coverage_mat) <- c(N_prioritisation_steps, N_age)

# Generating Vaccine Efficacy Over Time
vaccine_efficacy_infection_t[, ] <- interpolate(tt_vaccine_efficacy_infection, vaccine_efficacy_infection, "constant")
dim(vaccine_efficacy_infection_t) <- c(N_age, N_vaccine)
tt_vaccine_efficacy_infection[] <- user()
vaccine_efficacy_infection[, ,] <- user()
dim(tt_vaccine_efficacy_infection) <- user()
dim(vaccine_efficacy_infection) <- c(length(tt_vaccine_efficacy_infection), N_age, N_vaccine)

gamma_vaccine[] <- user() # Vector of vaccine progression parameters by vaccination status (First = 0 as handled separately as time-varying vaccination rate, Last = 0 as no progression from "previously vaccinated group)
dim(gamma_vaccine) <- N_vaccine

# Interpolation of vaccination rate over time
mv <- interpolate(tt_vaccine, max_vaccine, "constant")
tt_vaccine[] <- user()
max_vaccine[] <- user()
dim(tt_vaccine) <- user()
dim(max_vaccine) <- length(tt_vaccine)


# Track the proportion who have received vaccine in each age group

# Population size
pop_size[] <- sum(S[i,]) + sum(E1[i,]) + sum(E2[i,]) + sum(IMild[i,]) + sum(ICase1[i,]) + sum(ICase2[i,]) +
  sum(IMVGetLive1[i,]) + sum(IMVGetLive2[i,]) +
  sum(IMVGetDie1[i,]) + sum(IMVGetDie2[i,]) + sum(IMVNotGetLive1[i,]) + sum(IMVNotGetLive2[i,]) + sum(IMVNotGetDie1[i,]) + sum(IMVNotGetDie2[i,]) +
  sum(IOxGetLive1[i,]) + sum(IOxGetLive2[i,]) + sum(IOxGetDie1[i,]) + sum(IOxGetDie2[i,]) + sum(IOxNotGetLive1[i,]) + sum(IOxNotGetLive2[i,]) +
  sum(IOxNotGetDie1[i,]) + sum(IOxNotGetDie2[i,]) +
  sum(IRec1[i,]) + sum(IRec2[i,]) +
  sum(R1[i,]) + sum(R2[i,])
dim(pop_size) <- N_age
# Proportion who have received vaccine
pr[] <- 1 - ((sum(S[i,1]) + sum(E1[i,1]) + sum(E2[i,1]) + sum(IMild[i,1]) + sum(ICase1[i,1]) + sum(ICase2[i,1]) +
           sum(IMVGetLive1[i,1]) + sum(IMVGetLive2[i,1]) +
           sum(IMVGetDie1[i,1]) + sum(IMVGetDie2[i,1]) + sum(IMVNotGetLive1[i,1]) + sum(IMVNotGetLive2[i,1]) + sum(IMVNotGetDie1[i,1]) + sum(IMVNotGetDie2[i,1]) +
           sum(IOxGetLive1[i,1]) + sum(IOxGetLive2[i,1]) + sum(IOxGetDie1[i,1]) + sum(IOxGetDie2[i,1]) + sum(IOxNotGetLive1[i,1]) + sum(IOxNotGetLive2[i,1]) +
           sum(IOxNotGetDie1[i,1]) + sum(IOxNotGetDie2[i,1]) +
           sum(IRec1[i,1]) + sum(IRec2[i,1]) +
           sum(R1[i,1]) + sum(R2[i,1])) / (pop_size[i]))
dim(pr) <- N_age

# Isolate age groups below current target coverage which must be targeted
vaccination_target_mat[,] <- if(pr[j] < vaccine_coverage_mat[i,j]) 1 else 0
dim(vaccination_target_mat) <- c(N_prioritisation_steps, N_age)

vaccine_target_vec[] <- if(sum(vaccination_target_mat[i,]) == 0) 1 else 0
dim(vaccine_target_vec) <- N_prioritisation_steps
current_index <- min(sum(vaccine_target_vec) + 1, N_prioritisation_steps)

vaccination_target[] <- vaccination_target_mat[as.integer(current_index),i]
dim(vaccination_target) <- N_age

vr_temp[] <- S[i,1] * vaccination_target[i] + E1[i,1] * vaccination_target[i] + E2[i,1] * vaccination_target[i] + R1[i,1] * vaccination_target[i] + R2[i,1] * vaccination_target[i]
dim(vr_temp) <- N_age
# Catch so vaccination rate does not exceed 1 if the number of people available for vaccination < number of vaccines
vr_den <- if(sum(vr_temp) <= mv) mv else sum(vr_temp)
vr <- if(mv == 0) 0 else mv / vr_den  # Vaccination rate to achieve capacity given number in vaccine-eligible population
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

# Generating prob_hosp Over Time
prob_hosp_t[, ] <- interpolate(tt_vaccine_efficacy_disease, prob_hosp, "constant")
prob_hosp_t_mult[, ] <- prob_hosp_multiplier_t*prob_hosp_t[i,j]
dim(prob_hosp_t) <- c(N_age, N_vaccine)  # probability of requiring hospitalisation by age and vaccination status at time t
dim(prob_hosp_t_mult) <- c(N_age, N_vaccine)  # probability of requiring hospitalisation by age and vaccination status at time t
tt_vaccine_efficacy_disease[] <- user()
prob_hosp[, ,] <- user()
dim(tt_vaccine_efficacy_disease) <- user()
dim(prob_hosp) <- c(length(tt_vaccine_efficacy_disease), N_age, N_vaccine)

# Interpolation for prob_hosp_multiplier
prob_hosp_multiplier_t <- interpolate(tt_prob_hosp_multiplier, prob_hosp_multiplier, "constant")
tt_prob_hosp_multiplier[] <- user()
dim(tt_prob_hosp_multiplier) <- user()
dim(prob_hosp_multiplier) <- length(tt_prob_hosp_multiplier)
prob_hosp_multiplier[] <- user() # rate of progression through recovered compartment (loss of naturally acquired immunity)

# Probability of severe symptoms with multiplier
prob_severe[] <- user() # probability of severe disease (requiring mechanical ventilation) by age
dim(prob_severe) <- N_age

prob_severe_multiplier[] <- user() # flat modifer to severity
dim(prob_severe_multiplier) <- length(tt_prob_severe_multiplier)
tt_prob_severe_multiplier[] <- user()
dim(tt_prob_severe_multiplier) <- user()

#interpolate severity
prob_severe_multiplier_t <- interpolate(tt_prob_severe_multiplier, prob_severe_multiplier, "constant")
#calculate the new severity
prob_severe_multi[] <- prob_severe_multiplier_t*prob_severe[i]
dim(prob_severe_multi) <- N_age

prob_non_severe_death_treatment[] <- user() # probability of dying from non-severe disease (i.e. requiring oxygen but not mechanical ventilation) by age given you receive appropriate treatment (proxy here is whether a general hospital bed is available)
dim(prob_non_severe_death_treatment) <- N_age

prob_non_severe_death_no_treatment[] <- user() # probability of dying from non-severe disease (i.e. requiring oxygen but not mechanical ventilation) by age given you do NOT receive appropriate treatment (proxy here is whether a general hospital bed is available)
dim(prob_non_severe_death_no_treatment) <- N_age

prob_severe_death_treatment[] <- user() # probability of dying from severe disease (i.e. requiring mechanical ventilation) by age given you receive appropriate treatment (proxy here is whether an ICU bed is available)
dim(prob_severe_death_treatment) <- N_age

prob_severe_death_no_treatment[] <- user() # probability of dying from severe disease (i.e. requiring mechanical ventilation) by age given you do NOT receive appropriate treatment (proxy here is whether an ICU bed is available)
dim(prob_severe_death_no_treatment) <- N_age

rel_infectiousness[] <- user() # Relative infectiousness of age categories relative to maximum infectiousness age category
dim(rel_infectiousness) <- N_age

rel_infectiousness_vaccinated[,] <- user() # Relative infectiousness of vaccinated age categories relative to maximum infectiousness age category
dim(rel_infectiousness_vaccinated) <- c(N_age, N_vaccine)

# Infections Requiring Oxygen (a general Hosptial Bed)
hosp_occ <- sum(IOxGetLive1) + sum(IOxGetLive2) - gamma_get_ox_survive_t * sum(IOxGetLive2) + sum(IOxGetDie1) + sum(IOxGetDie2) - gamma_get_ox_die_t * sum(IOxGetDie2) + sum(IRec1) + sum(IRec2) - gamma_rec * sum(IRec2) # Summing number of infections in compartments that use general hospital beds
number_requiring_Ox[,] <- gamma_ICase * ICase2[i,j] * (1 - prob_severe_multi[i])
dim(number_requiring_Ox) <- c(N_age, N_vaccine)
total_number_requiring_ox <- sum(number_requiring_Ox)

p_oxygen <-  if ((total_number_requiring_ox <= (hosp_bed_capacity - hosp_occ)) || total_number_requiring_ox <= 0) 1 else (hosp_bed_capacity - hosp_occ) / total_number_requiring_ox

# Infections Requiring Mechanical Ventilation (an ICU Bed)
ICU_occ <- sum(IMVGetLive1) + sum(IMVGetLive2) - gamma_get_mv_survive_t * sum(IMVGetLive2) + sum(IMVGetDie1) + sum(IMVGetDie2) - gamma_get_mv_die_t * sum(IMVGetDie2) # Summing number of infections in compartments that use ICU beds
number_requiring_IMV[,] <- gamma_ICase * ICase2[i,j] * prob_severe_multi[i]
dim(number_requiring_IMV) <- c(N_age, N_vaccine)
total_number_requiring_IMV <- sum(number_requiring_IMV)

p_ventilation <- if (total_number_requiring_IMV <= (ICU_bed_capacity - ICU_occ) || total_number_requiring_IMV <= 0) 1 else (ICU_bed_capacity - ICU_occ) / total_number_requiring_IMV

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
temp_rel[,] <- (IMild[i,j] * rel_infectiousness_vaccinated[i,j]) + (ICase1[i,j] * rel_infectiousness_vaccinated[i,j]) + (ICase2[i,j] * rel_infectiousness_vaccinated[i,j])
temp[] <- sum(temp_rel[i,])
dim(temp_rel) <- c(N_age, N_vaccine)
dim(temp) <- c(N_age)

s_ij[,] <- m[i, j] * temp[j] * rel_infectiousness[j]
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

# Hospitalisations
deriv(hospitalisations_cumu[,]) <- number_requiring_IMV[i,j] * p_ventilation + number_requiring_Ox[i,j] * p_oxygen
dim(hospitalisations_cumu) <- c(N_age, N_vaccine)
initial(hospitalisations_cumu[,]) <- 0

# Deaths
output(deaths_cumu[,]) <- D[i,j]
dim(deaths_cumu) <- c(N_age, N_vaccine)

# Infections
deriv(infections_cumu[,]) <- (lambda[i] * vaccine_efficacy_infection_t[i,j] * S[i,j])
dim(infections_cumu) <- c(N_age, N_vaccine)
initial(infections_cumu[,]) <- 0

# Vaccinations
deriv(vaccines_cumu[]) <- (vr * vaccination_target[i] * S[i,1]) + (vr * vaccination_target[i] * E1[i,1]) + (vr * vaccination_target[i] * E2[i,1]) +
  (vr * vaccination_target[i] * R1[i,1]) + (vr * vaccination_target[i] * R2[i,1])
dim(vaccines_cumu) <- N_age
initial(vaccines_cumu[]) <- 0

# Unvaccinated
output(unvaccinated[]) <- sum(S[i,1]) + sum(E1[i,1]) + sum(E2[i,1]) + sum(IMild[i,1]) + sum(ICase1[i,1]) + sum(ICase2[i,1]) +
  sum(IMVGetLive1[i,1]) + sum(IMVGetLive2[i,1]) +
  sum(IMVGetDie1[i,1]) + sum(IMVGetDie2[i,1]) + sum(IMVNotGetLive1[i,1]) + sum(IMVNotGetLive2[i,1]) + sum(IMVNotGetDie1[i,1]) + sum(IMVNotGetDie2[i,1]) +
  sum(IOxGetLive1[i,1]) + sum(IOxGetLive2[i,1]) + sum(IOxGetDie1[i,1]) + sum(IOxGetDie2[i,1]) + sum(IOxNotGetLive1[i,1]) + sum(IOxNotGetLive2[i,1]) +
  sum(IOxNotGetDie1[i,1]) + sum(IOxNotGetDie2[i,1]) +
  sum(IRec1[i,1]) + sum(IRec2[i,1]) +
  sum(R1[i,1]) + sum(R2[i,1]) + sum(D[i,1])
dim(unvaccinated) <- N_age
# Vaccinated
output(vaccinated[]) <- sum(S[i,2:5]) + sum(E1[i,2:5]) + sum(E2[i,2:5]) + sum(IMild[i,2:5]) + sum(ICase1[i,2:5]) + sum(ICase2[i,2:5]) +
  sum(IMVGetLive1[i,2:5]) + sum(IMVGetLive2[i,2:5]) +
  sum(IMVGetDie1[i,2:5]) + sum(IMVGetDie2[i,2:5]) + sum(IMVNotGetLive1[i,2:5]) + sum(IMVNotGetLive2[i,2:5]) + sum(IMVNotGetDie1[i,2:5]) + sum(IMVNotGetDie2[i,2:5]) +
  sum(IOxGetLive1[i,2:5]) + sum(IOxGetLive2[i,2:5]) + sum(IOxGetDie1[i,2:5]) + sum(IOxGetDie2[i,2:5]) + sum(IOxNotGetLive1[i,2:5]) + sum(IOxNotGetLive2[i,2:5]) +
  sum(IOxNotGetDie1[i,2:5]) + sum(IOxNotGetDie2[i,2:5]) +
  sum(IRec1[i,2:5]) + sum(IRec2[i,2:5]) +
  sum(R1[i,2:5]) + sum(R2[i,2:5]) + sum(D[i,2:5])
dim(vaccinated) <- N_age
# Previously vaccinated
output(priorvaccinated[]) <- sum(S[i,6]) + sum(E1[i,6]) + sum(E2[i,6]) + sum(IMild[i,6]) + sum(ICase1[i,6]) + sum(ICase2[i,6]) +
  sum(IMVGetLive1[i,6]) + sum(IMVGetLive2[i,6]) +
  sum(IMVGetDie1[i,6]) + sum(IMVGetDie2[i,6]) + sum(IMVNotGetLive1[i,6]) + sum(IMVNotGetLive2[i,6]) + sum(IMVNotGetDie1[i,6]) + sum(IMVNotGetDie2[i,6]) +
  sum(IOxGetLive1[i,6]) + sum(IOxGetLive2[i,6]) + sum(IOxGetDie1[i,6]) + sum(IOxGetDie2[i,6]) + sum(IOxNotGetLive1[i,6]) + sum(IOxNotGetLive2[i,6]) +
  sum(IOxNotGetDie1[i,6]) + sum(IOxNotGetDie2[i,6]) +
  sum(IRec1[i,6]) + sum(IRec2[i,6]) +
  sum(R1[i,6]) + sum(R2[i,6]) + sum(D[i,6])
dim(priorvaccinated) <- N_age

output(N[]) <- pop_size[i] + sum(D[i,])
dim(N) <- N_age
################################################################################

