library(scales)
library(coda)
library(rjags)
library(BEST)
library(ggplot2)
library(tidyverse)
library(brms)
library(tidybayes)
library(NatParksPalettes)
library(jtools)
library(dplyr)


setwd("D:/Dissertation")

my_data <- read_csv('RdataImport_ENGR151.csv')


#Data Organization
Exp_Group <- my_data[my_data$Group == 'Exp',]
Con_Group <- my_data[my_data$Group == 'Control',]
Exp_cGPA <- Exp_Group[['Sem1gpa']]
Con_cGPA <- Con_Group[['Sem1gpa']]
Con_cGPA <- na.omit(Con_cGPA)

Exp_hGPA <- Exp_Group[['Hsgpa']]
Con_hGPA <- Con_Group[['Hsgpa']]

Exp_SAT <- Exp_Group[['SAT']]
Con_SAT <- Con_Group[['SAT']]
Exp_SAT <- na.omit(Exp_SAT)
Con_SAT <- na.omit(Con_SAT)

Exp_ACT <- Exp_Group[['ACT']]
Con_ACT <- Con_Group[['ACT']]
Exp_ACT <- na.omit(Exp_ACT)
Con_ACT <- na.omit(Con_ACT)

Pre_CAT <- Exp_Group[['CAT_pre']]
Post_CAT <- Exp_Group[['CAT_post']]
Delta_CAT <- Exp_Group[['CAT_delta']]
Pre_CAT <- na.omit(Pre_CAT)
Post_CAT <- na.omit(Post_CAT)
Delta_CAT <- na.omit(Delta_CAT)

muM_cGPA <- 2.82
muSD_cGPA <- 0.92

priors_cGPA <- list(muM = muM_cGPA, muSD = muSD_cGPA)

BESTcGPA <- BESTmcmc(Exp_cGPA, Con_cGPA, priors=priors_cGPA, parallel=FALSE)

plot(BESTcGPA)
plotAll(BESTcGPA)
#plot(1:100002, (BESTcGPA$mu1), type = "l")
#summary(BESTcGPA)
#print(BESTcGPA)

muM_hGPA <- 3.51
muSD_hGPA <- 0.46

priors_hGPA <- list(muM = muM_hGPA, muSD = muSD_hGPA)

BESThGPA <- BESTmcmc(Exp_hGPA, Con_hGPA, priors=priors_hGPA, parallel=FALSE)

plot(BESThGPA)
plotAll(BESThGPA)
#plot(1:100002, (BESThGPA$mu1), type = "l")
#summary(BESThGPA)
#print(BESThGPA)

muM_SAT <- 531
muSD_SAT <- 40

priors_SAT <- list(muM = muM_SAT, muSD = muSD_SAT)

BESTSAT <- BESTmcmc(Exp_SAT, Con_SAT, priors=priors_SAT, parallel=FALSE)

plot(BESTSAT)
plotAll(BESTSAT)
#plot(1:100002, (BESTSAT$mu1), type = "l")
#summary(BESTSAT)
#print(BESTSAT)

muM_ACT <- 24
muSD_ACT <- 1

priors_ACT <- list(muM = muM_ACT, muSD = muSD_ACT)

BESTACT <- BESTmcmc(Exp_ACT, Con_ACT, priors=priors_ACT, parallel=FALSE)

plot(BESTACT)
plotAll(BESTACT)
#plot(1:100002, (BESTACT$mu1), type = "l")
#summary(BESTACT)
#print(BESTACT)

muM_CAT <- 6.8
muSD_CAT <- 5

priors_CAT <- list(muM = muM_CAT, muSD = muSD_CAT)

BESTCAT <- BESTmcmc(Post_CAT, Pre_CAT, priors=priors_CAT, parallel=FALSE)
BESTCATD <- BESTmcmc(Delta_CAT, priors=priors_CAT, parallel=FALSE)

plot(BESTCAT)
plotAll(BESTCAT)

plot(BESTCATD)
plotAll(BESTCATD)

#___________________________________________________________________
#Proportion Testing

library(tidyverse)    # ggplot, dplyr, and friends
library(gt)           # Fancy tables
library(glue)         # Easier string interpolation
library(scales)       # Nicer labeling functions
library(ggmosaic)     # Mosaic plots with ggplot
library(ggpattern)    # Pattern fills in ggplot
library(patchwork)    # Combine plots nicely
library(parameters)   # Extract model parameters as data frames
library(cmdstanr)     # Run Stan code from R
library(brms)         # Nice frontend for Stan
library(tidybayes)    # Manipulate MCMC chains in a tidy way
library(likert)       # Contains the pisaitems data
library(jtools)

options(brms.backend = "cmdstanr")

clrs_saguaro <- NatParksPalettes::natparks.pals("Saguaro")
clr_Exp <- clrs_saguaro[1]
clr_Con <- clrs_saguaro[6]
clr_Diff <- clrs_saguaro[4]

CHAINS <- 4
ITER <- 3000
WARMUP <- 1500
BAYES_SEED <- 1231

A_counts <- my_data %>%
  group_by(Group, Alg_A) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

A_values <- A_counts %>%
  filter(Alg_A == 1)

A_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = A_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  #need to up the warm up and iterations because there seemed to be a
  #convergence issue
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


Ap1 <- A_model %>% 
  epred_draws(newdata = A_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students with an A in College Algebra",
       y = NULL) +
  theme_nice()

Ap2 <- A_model %>% 
  epred_draws(newdata = A_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "pp." isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the pointrange doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(Ap1 / plot_spacer() / Ap2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

A_delta <- A_model %>% 
  epred_draws(newdata = A_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)
#___________________

DFW_counts <- my_data %>%
  group_by(Group, Alg_DFW) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

DFW_values <- DFW_counts %>%
  filter(Alg_DFW == 1)

DFW_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = DFW_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


DFWp1 <- DFW_model %>% 
  epred_draws(newdata = DFW_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015),
                     limits = c(0,0.5)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students with a DFW in College Algebra",
       y = NULL) +
  theme_nice()

DFWp2 <- DFW_model %>% 
  epred_draws(newdata = DFW_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(DFWp1 / plot_spacer() / DFWp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

DFW_delta <- DFW_model %>% 
  epred_draws(newdata = DFW_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

#___________________

ENGRpersist_counts <- my_data %>%
  group_by(Group, ENGR_persist) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

EP_values <- ENGRpersist_counts %>%
  filter(ENGR_persist == 1)

EP_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = EP_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


EPp1 <- EP_model %>% 
  epred_draws(newdata = EP_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students who persisted in Engineering",
       y = NULL) +
  theme_nice()

EPp2 <- EP_model %>% 
  epred_draws(newdata = EP_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(EPp1 / plot_spacer() / EPp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

EP_delta <- EP_model %>% 
  epred_draws(newdata = EP_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

#___________________

Ipersist_counts <- my_data %>%
  group_by(Group, INST_persist) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

IP_values <- Ipersist_counts %>%
  filter(INST_persist == 1)

IP_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = IP_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


IPp1 <- IP_model %>% 
  epred_draws(newdata = IP_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students who persisted in the Institution",
       y = NULL) +
  theme_nice()

IPp2 <- IP_model %>% 
  epred_draws(newdata = IP_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(IPp1 / plot_spacer() / IPp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

IP_delta <- IP_model %>% 
  epred_draws(newdata = IP_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

#___________________

OT2_counts <- my_data %>%
  group_by(Group, On_Track_2) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

OT2_values <- OT2_counts %>%
  filter(On_Track_2 == 1)

OT2_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = OT2_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


OT2p1 <- OT2_model %>% 
  epred_draws(newdata = OT2_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of Students On Track: Sem 2",
       y = NULL) +
  theme_nice()

OT2p2 <- OT2_model %>% 
  epred_draws(newdata = OT2_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(OT2p1 / plot_spacer() / OT2p2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

OT2_delta <- OT2_model %>% 
  epred_draws(newdata = OT2_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

#___________________

OT3_counts <- my_data %>%
  group_by(Group, On_Track_3) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

OT3_values <- OT3_counts %>%
  filter(On_Track_3 == 1)

OT3_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = OT3_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


OT3p1 <- OT3_model %>% 
  epred_draws(newdata = OT3_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of Students On Track: Sem 3",
       y = NULL) +
  theme_nice()

OT3p2 <- OT3_model %>% 
  epred_draws(newdata = OT3_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(OT3p1 / plot_spacer() / OT3p2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

OT3_delta <- OT3_model %>% 
  epred_draws(newdata = OT3_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

#___________________

OT4_counts <- my_data %>%
  group_by(Group, On_Track_4) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

OT4_values <- OT4_counts %>%
  filter(On_Track_4 == 1)

OT4_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = OT4_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


OT4p1 <- OT4_model %>% 
  epred_draws(newdata = OT4_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of Students On Track: Sem 4",
       y = NULL) +
  theme_nice()

OT4p2 <- OT4_model %>% 
  epred_draws(newdata = OT4_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(OT4p1 / plot_spacer() / OT4p2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

OT4_delta <- OT4_model %>% 
  epred_draws(newdata = OT4_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

#___________________

NM4_counts <- my_data %>%
  group_by(Group, No_Math_4) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

NM4_values <- NM4_counts %>%
  filter(No_Math_4 == 1)

NM4_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = NM4_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


NM4p1 <- NM4_model %>% 
  epred_draws(newdata = NM4_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of Students to Leave Track: Sem 4",
       y = NULL) +
  theme_nice()

NM4p2 <- NM4_model %>% 
  epred_draws(newdata = NM4_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "%" isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the point range doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(NM4p1 / plot_spacer() / NM4p2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

NM4_delta <- NM4_model %>% 
  epred_draws(newdata = NM4_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)
