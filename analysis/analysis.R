# ============================================================================ #
# ** PRELIMINARIES ** ####
# ============================================================================ #

# load packages ...

if(!require('pacman'))install.packages('pacman')

pacman::p_load(tidyverse, 
               lmerTest, 
               beepr, 
               tictoc, 
               ggh4x, 
               emmeans, 
               performance, 
               ordinal, 
               kableExtra, 
               patchwork, 
               gridExtra, 
               ggtext)

# load data ...

df <- read_csv('data/data.csv')

# define functions ...

my.lcl <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$conf.int[1])
}

my.ucl <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$conf.int[2])
}

my_table <- function(df) {
  if ("Estimate" %in% colnames(df)) {
    df <- df %>% rename("Estimate~($\\hat{\\beta}$)" = "Estimate")
  }
  if ("estimate" %in% colnames(df)) {
    df <- df %>% rename("Estimate~($\\hat{\\beta}$)" = "estimate")
  }
  if ("Std. Error" %in% colnames(df)) {
    df <- df %>% rename("$SE$" = "Std. Error")
  }
  if ("t value" %in% colnames(df)) {
    df <- df %>%
      rename("$t$" = "t value")
  }
  if ("t.ratio" %in% colnames(df)) {
    df <- df %>%
      rename("$t$" = "t.ratio")
  }
  if ("z value" %in% colnames(df)) {
    df <- df %>% rename("$z$" = "z value")
  }
  if ("z.ratio" %in% colnames(df)) {
    df <- df %>%
      rename("$z$" = "z.ratio")
  }
  if ("Pr(>|z|)" %in% colnames(df)) {
    df <- df %>%
      rename("$p$" = "Pr(>|z|)") %>%
      mutate(`$p$` = scales::pvalue(`$p$`))
  }
  if ("Pr(>|t|)" %in% colnames(df)) {
    df <- df %>%
      rename("$p$" = "Pr(>|t|)") %>%
      mutate(`$p$` = scales::pvalue(`$p$`))
  }
  if ("p.value" %in% colnames(df)) {
    df <- df %>%
      rename("$p$" = "p.value") %>%
      mutate(`$p$` = scales::pvalue(`$p$`)) 
  }
  if ("contrast" %in% colnames(df)) {
    df <- df %>%
      rename("Contrast" = "contrast")
  }
  if ("SE" %in% colnames(df)) {
    df <- df %>%
      rename("$SE$" = "SE")
  }
  if ("condition" %in% colnames(df)) {
    df <- df %>%
      rename("Condition" = "condition")
  }
  if ("group" %in% colnames(df)) {
    df <- df %>%
      rename("Group" = "group")
  }
  if ("environment" %in% colnames(df)) {
    df <- df %>%
      rename("Environment" = "environment")
  }
  if ("dependency" %in% colnames(df)) {
    df <- df %>%
      rename("Dependency" = "dependency")
  }
  if ("pairing" %in% colnames(df)) {
    df <- df %>%
      rename("Pairing" = "pairing")
  }
  if ("mean" %in% colnames(df)) {
    df <- df %>%
      rename("Mean" = "mean")
  }
  if ("lcl" %in% colnames(df)) {
    df <- df %>%
      rename("LCL" = "lcl")
  }
  if ("ucl" %in% colnames(df)) {
    df <- df %>%
      rename("UCL" = "ucl")
  }
  if ("prop" %in% colnames(df)) {
    df <- df %>%
      rename("\\%" = "prop")
  }
  if ("proficiency.trend" %in% colnames(df)) {
    df <- df %>%
      rename("Proficiency Trend" = "proficiency.trend")
  }
  if ("region2" %in% colnames(df)) {
    df <- df %>%
      mutate(region2 = as.character(region2)) %>%
      rename("Region" = "region2")
  }
  if ("mean_rt" %in% colnames(df)) {
    df <- df %>%
      rename("Mean" = "mean_rt")
  }
  df <- df %>%
    #mutate(Count = round(Count)) %>%
    #mutate(Count = as.character(Count)) %>%
    mutate_if(is.numeric, ~ format(round(., 2), nsmall = 2)) %>%
    mutate_all(as.character) %>%
    mutate_if(is.character, str_replace_all, 
              c('-' = '−', 
                '<0.001' = '< .001',
                '>0.999' = '> .999',
                'group2' = 'Group[KLE]', 
                'group3' = 'Group[MLE]',
                'dependency2' = 'Dependency[RP]',
                'environment2' = 'Environment[Long]',
                'environment3' = 'Environment[Island]',
                'cond1' = 'Short',
                'cond2' = 'Long',
                'cond3' = 'Island',
                'short' = 'Short', 
                'long' = 'Long', 
                'island' = 'Island',
                'gap' = 'Gap',
                'resumption' = 'Resumption',
                'pronoun' = 'RP',
                'other' = 'Other',
                'proficiency' = 'Proficiency',
                'grammatical' = 'Grammatical',
                'unGrammatical' = 'Ungrammatical',
                'pairing2' = 'Pairing[KE]', 
                'pairing3' = 'Pairing[ME]',
                'pairing4' = 'Pairing[KK]',
                'pairing5' = 'Pairing[MM]',
                'set1' = 'F1',
                'set2' = 'F2',
                'set3' = 'F3',
                'set4' = 'F4',
                'set5' = 'F5',
                'set6' = 'F6')) %>%
    mutate_if(is.character, str_replace_all,
              c(':Environment\\[Long\\]:' = ':Envir.\\[Long\\]:',
                ':Environment\\[Island\\]:' = ':Envir.\\[Island\\]:'))
    if ("LCL" %in% colnames(df)) {
      df <- df %>%
        mutate(LCL = paste0("[", LCL, ", ", UCL, "]")) %>%
        rename(`95\\% $CI$` = LCL) %>%
        select(-UCL)
    }
  return(df)
}

csv_table <- function(df) {
  if ("Estimate~($\\hat{\\beta}$)" %in% colnames(df)) {
    df <- df %>% rename("Estimate" = "Estimate~($\\hat{\\beta}$)")
  }
  if ("\\%" %in% colnames(df)) {
    df <- df %>%
      rename("Percent" = "\\%")
  }
  if ("%" %in% colnames(df)) {
    df <- df %>%
      rename("Percent" = "%")
  }
  if ("Participants (\\%)" %in% colnames(df)) {
    df <- df %>%
      rename("Percent" = "Participants (\\%)")
  }
  if ("Accuracy (\\%)" %in% colnames(df)) {
    df <- df %>%
      rename("Percent" = "Accuracy (\\%)")
  }
  if ("95\\% $CI$" %in% colnames(df)) {
    df <- df %>%
      rename("$CI$" = "95\\% $CI$")
  }
  return(df)
}

# ============================================================================ #
# ** COUNT PARTICIPANTS BY TASK ** ####
# ============================================================================ #

# plot number of participants in each task  ...

plot <- df %>%
  group_by(study, group, task, participant) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study, group, task) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(task = as.factor(task), group = as.factor(group))

ggplot(data = plot, aes(x = task, y = n, group = group, col = group, shape = group, label = n)) +
  geom_line() +
  geom_point(lwd = 3, alpha = .5) +
  geom_text(size = 3, col = 'black', nudge_y = 1.5) +
  theme_classic() +
  scale_x_discrete(name='Task', 
                   limits = c('language_survey', 'ept', 'sprt', 'english_ajt', 'korean_ajt', 'mandarin_ajt', 'ctest', 'exit_survey'), 
                   labels = c('LBQ', 'EPT', 'SPRT', 'EAJT', 'KAJT', 'MAJT', 'C-test', 'Exit')) +
  scale_y_continuous(name='Participants', limits=c(50, 100)) +
  scale_color_manual(name = 'Group', values = c('#FFC107', '#1E88E5', '#D81B60')) +
  scale_shape_manual(name = 'Group', values = c(16, 15, 17)) +
  theme(text = element_text(size = 12),
        legend.position = 'bottom') +
  facet_wrap(~study)

ggsave("plots/participant_counts.png", width=6.5, height=4.5, dpi=600)

remove(plot)

# ============================================================================ #
# ** LANGUAGE BACKGROUND SURVEY ** ####
# ============================================================================ #

# prep dataframe ...

survey <- df %>%
  filter(task == 'language_survey') %>%
  arrange(study, group, participant) %>%
  select_if(function(x){!all(is.na(x))})

remove(survey)

# ============================================================================ #
# ** ENGLISH PROFICIENCY C-TEST ** ####
# ============================================================================ #

# prep dataframe ...

ct <- df %>% filter(task == 'ctest') %>%
  arrange(study, group, participant) %>%
  select_if(function(x){!all(is.na(x))})

# ---------------------------------------------------------------------------- #
# plot: accuracy ####
# ---------------------------------------------------------------------------- #

# summarise ...

plot <- df %>%
  group_by(study, group, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  filter(is.na(proficiency) == FALSE) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant)),
         mean = mean(proficiency),
         sd = sd(proficiency),
         min = min(proficiency),
         max = max(proficiency)) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')'))

# check participants with low accuracy ...

check <- plot %>%
  mutate(low = case_when(proficiency < 25 ~ TRUE, TRUE ~ FALSE)) %>%
  group_by(study, group) %>%
  summarise(prop_low = mean(low),
            n_low = length(which(low == TRUE)),
            n_group = length(unique(participant)))

# define plot data ...

p <- ggplot(data=plot)
p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# define plot styling ...

s <- list(
  aes(x=study, y=proficiency, col=group, fill=group),
  geom_richtext(y = 2, 
                col = "black",
                hjust = .5,
                size = 2.95,
                fill = NA, 
                label.color = NA,
                aes(label = paste0("Mean = ", 
                                   format(round(mean, 2), nsmall = 2),
                                   ", *SD* = ",
                                   format(round(sd, 2), nsmall = 2),
                                   ", ",
                                   min,
                                   "--",
                                   max
                                   ))),
  geom_hline(yintercept = 25),
  geom_violin(),
  geom_boxplot(width = .1, fill = 'white'),
  theme_classic(),
  xlab('Group'),
  scale_y_continuous(name='Score',
                     limits=c(0, 50),
                     breaks = c(0, 25, 50)),
  scale_color_manual(name = 'Group', 
                    values = c('coral', '#ffb000', '#40B0A6')),
  scale_fill_manual(name = 'Group', 
                    values = c('#ffccb9', '#ffdf99', '#b0e2de')),
  theme(text = element_text(size = 12),
        legend.position = 'hide',
        plot.margin = margin(0, 0, 0, 0),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_markdown()
        ),
  facet_grid(study~group_label)
)

# print and save plots ...

p1 + s
ggsave("plots/orc_ctest.png", width=6.5, height=1.5, dpi=600)

p2 + s
ggsave("plots/src_ctest.png", width=6.5, height=1.5, dpi=600)

# check dispersion ...

check <- plot %>%
  group_by(study, group) %>%
  summarise(mean = mean(proficiency, na.rm = TRUE),
            sd = sd(proficiency, na.rm = TRUE),
            iqr = IQR(proficiency, na.rm = TRUE)) %>%
  ungroup()
  
remove(check, plot, s, p, p1, p2)

# ---------------------------------------------------------------------------- #
# model: accuracy ####
# ---------------------------------------------------------------------------- #

# filter for analysis ...

md <- ct %>%
  select(study, group, participant, item, accuracy) %>%
  mutate_at(c('study', 'group', 'participant', 'item'), factor)

# set contrasts for deviation coding ...

contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

# check contrasts ...

contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model <- glmer(accuracy ~ group + 
                  (1 | participant) + 
                  (1 + group | item), 
                data = filter(md, study == 'ORC'), 
                family = binomial, 
                control = glmerControl(optimizer = "bobyqa"))
summary(model)
toc()
beep()

# 35 sec elapsed
# no warnings

# save ...

model %>% write_rds('models/orc_ctest.rds')
model <- read_rds('models/orc_ctest.rds')

# performance ...

check_convergence(model) # okay
check_singularity(model) # okay
MuMIn::r.squaredGLMM(model) # theoretical: R2m = 0.1595369, R2c = 0.6370303; delta: R2m = 0.1430333 , R2c = 0.5711313

# model summary ...

summary(model)

# Parameter   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.3629     0.2464   5.531 3.18e-08 ***
# group2       -2.4521     0.2881  -8.511  < 2e-16 ***
# group3       -2.5051     0.2943  -8.513  < 2e-16 ***

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/orc_ctest_accuracy_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ctest_accuracy_model.csv')

# post-hoc testing ...

pairwise <- model %>%
  emmeans(~ group) %>%
  contrast('revpairwise', simple = 'group', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast  estimate    SE  df z.ratio p.value
# KLE - ENS   -2.452 0.288 Inf  -8.511  <.0001 ***
# MLE - ENS   -2.505 0.294 Inf  -8.513  <.0001 ***
# MLE - KLE   -0.053 0.187 Inf  -0.283  0.9557
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: mvt method for 3 tests 

# table ...

dt <- pairwise %>%
  my_table() %>%
  select(-df)

dt %>% write_rds('tables/orc_ctest_accuracy_posthoc.rds')

csv <- dt %>% 
  csv_table() %>% 
  write_csv('tables/orc_ctest_accuracy_posthoc.csv')

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model <- glmer(accuracy ~ group + 
                  (1 | participant) + 
                  (1 + group | item), 
                data = filter(md, study == 'SRC'), 
                family = binomial, 
                control = glmerControl(optimizer = "bobyqa"))
summary(model)
toc()
beep()

# 23 sec elapsed
# no warnings

model %>% write_rds('models/src_ctest.rds')
model <- read_rds('models/src_ctest.rds')

# model summary ...

summary(model)

# Parameter   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.5141     0.2463   6.147 7.91e-10 ***
# group2       -2.1268     0.3072  -6.924 4.39e-12 ***
# group3       -2.5014     0.3253  -7.689 1.49e-14 ***

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_ctest_accuracy_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/src_ctest_accuracy_model.csv')
dt <- read_rds('tables/src_ctest_accuracy_model.rds')

# post-hoc testing ...

pairwise <- model %>%
  emmeans(~ group) %>%
  contrast('revpairwise', simple = 'group', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast  estimate    SE  df z.ratio p.value
# KLE - ENS   -2.127 0.307 Inf  -6.924  <.0001 ***
# MLE - ENS   -2.501 0.325 Inf  -7.689  <.0001 ***
# MLE - KLE   -0.375 0.204 Inf  -1.836  0.1528
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: mvt method for 3 tests

# table ...

dt <- pairwise %>%
  select(contrast, estimate, SE, z.ratio, p.value) %>%
  my_table()

dt %>% write_rds('tables/src_ctest_accuracy_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ctest_accuracy_posthoc.csv')
dt <- read_rds('tables/src_ctest_accuracy_posthoc.rds')

remove(dt, pairwise, csv, model, md)

# ============================================================================ #
# ** ELICITED PRODUCTION TASK ** ####
# ============================================================================ #

# prep dataframe ...

ep <- df %>%
  filter(task == 'ept') %>%
  arrange(study, group, participant) %>%
  select_if(function(x){!all(is.na(x))})

check <- ep %>%
  group_by(study, group) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

# remove participants who gave only nontarget responses ...

list <- ep %>%
  filter(type != 'nontarget') %>%
  group_by(study, group, participant) %>%
  summarise() %>%
  ungroup()

check <- ep %>%
  group_by(study, group, participant) %>%
  summarise() %>%
  ungroup() %>%
  filter(!participant %in% c(list$participant))

ep <- ep %>%
  filter(participant %in% list$participant)

remove(check, list)

# ---------------------------------------------------------------------------- #
# plot: token critical (gap, resumption, and other responses) ####
# ---------------------------------------------------------------------------- #

# summarise for plotting ...

plot <- ep %>%
  filter(is.na(condition) == FALSE) %>%
  mutate(type = as.character(type)) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study, group, participant, condition, type, group_label) %>%
  count(type) %>%
  ungroup() %>%
  complete(type, nesting(study, group, participant, condition, group_label), fill = list(n = 0)) %>%
  select(study, group, participant, condition, type, n, group_label) %>%
  arrange(study, group, participant, condition) %>%
  group_by(study, group, condition, type, group_label) %>%
  summarise(
    mean = mean(n, na.rm = TRUE),
    lcl = my.lcl(n, conf.level = .95),
    ucl = my.ucl(n, conf.level = .95)
    ) %>%
  mutate(lcl = case_when(lcl == 'NaN' ~ 0, TRUE ~ lcl),
         ucl = case_when(ucl == 'NaN' ~ 0, TRUE ~ ucl),
         prop = mean * 20
         ) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = c('gap', 'resumption', 'other', 'nontarget'))) %>%
  filter(type %in% c('gap', 'resumption', 'other'))

# define data for plots ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# styling ...

s <- list(
  aes(x=condition, y=mean, group=type, col=type, shape=type, label=paste0(round(prop))),
  geom_hline(yintercept = 2.5),
  geom_line(aes(linetype = type, group = type), position = position_dodge(width = .6)),
  # scale_colour_manual(name='Dependency', 
  #                     values=c('#9ab6ff', '#ffcc59'),
  #                     breaks = c('gap', 'resumption'),
  #                     labels = c('Gap', 'Resumption')),
  # new_scale('colour'),
  geom_errorbar(aes(ymin = lcl, ymax = ucl, col = type), 
                linetype = 1,
                #width = .5,
                width = .2,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(aes(col = type), lwd = 1, position = position_dodge(width = .6), shape = 21, 
             fill = 'white', show.legend = FALSE),
  theme_classic(),
  #geom_text(size = 2.5, col = "black", hjust = 0, vjust = -.5, position = position_dodge(width = .6)),
  scale_y_continuous(name='Tokens', 
                     limits = c(-.05, 5),
                     breaks = c(0, 2.5, 5)),
  scale_x_discrete(name='Environment', 
                   labels = c('Short', 'Long', 'Island')
  ),
  scale_linetype_manual(name = 'Dependency:',
                        values = c('41', 'solid', '7222'),
                        breaks = c('gap', 'resumption', 'other'),
                        labels = c('Gapping', 'Resumption', 'Other')),
  scale_colour_manual(name='Dependency:', 
                      values=c('#648fff', '#ffb000', 'coral'),
                      breaks = c('gap', 'resumption', 'other'),
                      labels = c('Gapping', 'Resumption', 'Other')),
  theme(text = element_text(size = 12),
        plot.margin = margin(2, 2, 0, 2),
        legend.position = 'right', # dissertation
        #legend.position = 'bottom', # bucld
        legend.margin = margin(t = -.4, unit = 'cm'),
        #legend.key.width = unit(.8, 'cm'),
        plot.title = element_text(size = 12, hjust = .5),
        strip.text = element_markdown()),
  facet_manual(~group_label,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ept_token.png', width=6.5, height=2, dpi=600)

p2 + s
ggsave('plots/src_ept_token.png', width=6.5, height=2, dpi=600)

remove(plot, s, p1, p2)

# ---------------------------------------------------------------------------- #
# plot: token filler ####
# ---------------------------------------------------------------------------- #

# summarise for plotting ...

plot <- ep %>%
  filter(is.na(condition) == TRUE) %>%
  mutate(type = as.character(type)) %>%
  mutate(set = case_when(item %in% c('item16', 'item17', 'item18', 'item19', 'item20') ~ 'set1',
                         item %in% c('item21', 'item22', 'item23', 'item24', 'item25') ~ 'set2',
                         item %in% c('item26', 'item27', 'item28', 'item29', 'item30') ~ 'set3')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  ungroup() %>%
  group_by(study, group, participant, set, type, group_label) %>%
  count(type) %>%
  ungroup() %>%
  complete(type, nesting(study, group, participant, set, group_label), fill = list(n = 0)) %>%
  select(study, group, participant, set, type, n, group_label) %>%
  arrange(study, group, participant, set) %>%
  group_by(study, group, set, type, group_label) %>%
  summarise(
    mean = mean(n, na.rm = TRUE),
    lcl = my.lcl(n, conf.level = .95),
    ucl = my.ucl(n, conf.level = .95)
  ) %>%
  mutate(lcl = case_when(lcl == 'NaN' ~ 0, TRUE ~ lcl),
         ucl = case_when(ucl == 'NaN' ~ 0, TRUE ~ ucl),
         prop = mean * 20
  ) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = c('gap', 'other', 'nontarget', 'resumption'))) %>%
  filter(!type %in% c('nontarget')) %>%
  mutate(type = fct_relevel(type, "gap", "resumption", "other"))

# define data for plots ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# styling ...

s <- list(
  aes(x=set, y=mean, group=type, col=type, shape=type, label=paste0(round(prop), '%')),
  geom_hline(yintercept = 2.5),
  geom_line(aes(linetype = type, group = type), position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(aes(col = type), lwd = 1, position = position_dodge(width = .6), shape = 21,
             fill = 'white', show.legend = FALSE),
  theme_classic(),
  scale_y_continuous(name='Tokens', limits = c(-.05, 5),
                     breaks = c(0, 2.5, 5)),
  scale_x_discrete(name = 'Filler Type', labels = c('F1', 'F2', 'F3')),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000', 'coral'),
                      breaks = c('gap', 'resumption', 'other'),
                      labels = c('Gapping', 'Resumption', 'Other')),
  scale_linetype_manual(name='Dependency', 
                        values=c('41', 'solid', '7222'),
                        breaks = c('gap', 'resumption', 'other'),
                        labels = c('Gapping', 'Resumption', 'Other')),
  theme(text = element_text(size = 12),
        legend.position = 'right',
        legend.margin = margin(t = -.4, unit = 'cm'),
        strip.text = element_markdown()),
  facet_manual(~group_label,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ept_token_filler.png', width=6.5, height=2, dpi=600)

p2 + s
ggsave('plots/src_ept_token_filler.png', width=6.5, height=2, dpi=600)

# table ...

dt <- plot %>%
  select(-group_label) %>%
  arrange(study, group, set, type) %>%
  mutate(type = as.character(type)) %>%
  rename('Type' = 'set',
         'Dependency' = 'type', 
         'Tokens' = 'mean') %>%
  mutate_if(is.character, str_replace_all, 
            c('set1' = 'F1',
              'set2' = 'F2',
              'set3' = 'F3',
              'gap' = 'Gapping')) %>%
  my_table()

dt %>% write_rds('tables/ept_token_filler.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/ept_token_filler.csv')
dt <- read_rds('tables/ept_token_filler.rds')

remove(p1, p2, s, plot, csv)

# ---------------------------------------------------------------------------- #
# model: resumption ####
# ---------------------------------------------------------------------------- #

# prep data for orc modeling ...

md1 <- ep %>%
  mutate(resumption = case_when(type == 'resumption' ~ TRUE,
                                TRUE ~ FALSE)) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(group = fct_drop(group)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  select(study, group, participant, item, environment, resumption, type) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'environment', 'type'), factor) %>%
  filter(is.na(environment) == FALSE) %>%
  mutate(environment = fct_drop(environment))

summary(md1)

contrasts(md1$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md1$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md1$environment)
contrasts(md1$group)

# prep data for src modeling ...

md2 <- md1 %>%
  filter(environment %in% c('long', 'island')) %>%
  mutate(environment = fct_drop(environment))

summary(md2)

contrasts(md2$environment) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md2$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md2$environment)
contrasts(md2$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model <- glmer(resumption ~ environment * group +
                   (1 | participant) + 
                   (1 + environment + group | item), 
                 data = filter(md1, study == 'ORC'), family = binomial, 
                 control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
summary(model)
toc()
beep()

# ... 196.86 sec elapsed

# check performance ...

check_convergence(model) # okay
check_singularity(model, tolerance = 1e-9) # warning

# save ...

model %>% write_rds('models/orc_ept_resumption_model.rds')
model <- read_rds('models/orc_ept_resumption_model.rds')

# model summary ...

summary(model)

# Parameter            Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -5.40986    0.41999 -12.881  < 2e-16 ***
# environment2         1.75872    0.40678   4.324 1.54e-05 ***
# environment3         3.98097    0.39969   9.960  < 2e-16 ***
# group2               0.02127    0.69208   0.031  0.97548    
# group3              -0.20318    0.71451  -0.284  0.77613    
# environment2:group2  0.55035    0.87442   0.629  0.52909    
# environment3:group2 -2.78794    0.86416  -3.226  0.00125 ** 
# environment2:group3  0.78253    0.94940   0.824  0.40981    
# environment3:group3 -3.09904    0.94242  -3.288  0.00101 ** 

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% kbl(format = 'markdown')

# |Parameter                      | Estimate|   SE|      z|     p|
# |:------------------------------|--------:|----:|------:|-----:|
# |(Intercept)                    |    -5.41| 0.42| -12.88| 0.000| ***
# |Environment[Long]              |     1.76| 0.41|   4.32| 0.000| ***
# |Environment[Island]            |     3.98| 0.40|   9.96| 0.000| ***
# |Group[KLE]                     |     0.02| 0.69|   0.03| 0.975|
# |Group[MLE]                     |    -0.20| 0.71|  -0.28| 0.776|
# |Environment[Long]:Group[KLE]   |     0.55| 0.87|   0.63| 0.529|
# |Environment[Island]:Group[KLE] |    -2.79| 0.86|  -3.23| 0.001| **
# |Environment[Long]:Group[MLE]   |     0.78| 0.95|   0.82| 0.410|
# |Environment[Island]:Group[MLE] |    -3.10| 0.94|  -3.29| 0.001| **

dt %>% write_rds('tables/orc_ept_resumption_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ept_resumption_model.csv')

# joint tests ...

joint_tests(model)

# model term        df1 df2 F.ratio p.value
# environment         2 Inf  56.867  <.0001
# group               2 Inf   0.060  0.9421
# environment:group   4 Inf  11.664  <.0001

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ environment * group) %>%
  contrast('revpairwise', simple = 'environment', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast       group estimate    SE  df z.ratio p.value
# long - short   ENS      1.314 0.720 Inf   1.827  0.3981
# island - short ENS      5.943 0.717 Inf   8.289  <.0001
# island - long  ENS      4.629 0.510 Inf   9.080  <.0001
# long - short   KLE      1.865 0.558 Inf   3.342  0.0069
# island - short KLE      3.155 0.553 Inf   5.705  <.0001
# island - long  KLE      1.291 0.420 Inf   3.074  0.0173
# long - short   MLE      2.097 0.664 Inf   3.159  0.0131
# island - short MLE      2.844 0.647 Inf   4.394  0.0001
# island - long  MLE      0.747 0.444 Inf   1.683  0.4984

# table ...

dt <- pairwise %>%
  select(group, contrast, estimate, SE, z.ratio, p.value) %>%
  my_table()

dt %>% write_rds('tables/orc_ept_resumption_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ept_resumption_posthoc.csv')

# plot ...

plot(emmeans(model, ~ environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

remove(csv, model, dt)

# ============================================================================ #
# ... src model (final) ####
# ============================================================================ #

# src resumption rates with dummy participants

# load SRC elicited production task data ...

# add dummy participants (with only resumption responses) to ENS and KLE groups ...

dummy <- read_csv('data/ept_dummy.csv') %>% 
  filter(group %in% c('ENS', 'KLE'))

md <- bind_rows(md2, dummy) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'environment', 'type'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(study = fct_drop(study)) %>%
  mutate(environment = fct_drop(environment))

# add proficiency scores ...

prof <- read_rds("data/proficiency.rds")

md <- md %>%
  left_join(prof)

md <- md %>%
  mutate_at(c('study', 'group', 'participant'), factor)

# inspect dataframe ...

summary(md)

# set contrasts ...

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$environment)
contrasts(md$group)

# fit model ...

tic()
model <- glmer(resumption ~ environment * group +
                  (1 + environment | participant) + 
                  (1 | item), 
                data = md, family = binomial, 
                control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
summary(model)
toc()
beep()

# fixed-effects output:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -5.78620    1.27232  -4.548 5.42e-06 ***
# environment2         9.34694    3.66671   2.549 0.010799 *  
# environment3        14.09486    3.65409   3.857 0.000115 ***
# group2              -0.07404    1.82196  -0.041 0.967587    
# group3               1.94617    1.53277   1.270 0.204190    
# environment2:group2  6.60606    4.79703   1.377 0.168477    
# environment3:group2 -2.57811    4.82118  -0.535 0.592824    
# environment2:group3  7.03231    3.86099   1.821 0.068550 .  
# environment3:group3 -4.22413    3.83685  -1.101 0.270925 

# save ...

model %>% write_rds('models/src_ept_resumption_model.rds')
model <- read_rds('models/src_ept_resumption_model.rds')

# table ...

summary(model)$coefficients

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  # mutate_if(is.character, str_replace_all, 
  #           c('group2' = 'group3')) %>%
  my_table()

dt %>% write_rds('tables/src_ept_resumption_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ept_resumption_model.csv')
dt <- read_rds('tables/src_ept_resumption_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                      |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:------------------------------|:------------------------|:----|:-----|:------|
# |(Intercept)                    |−5.79                    |1.27 |−4.55 |< .001 |
# |Environment[Island]            |9.35                     |3.67 |2.55  |0.011  |
# |Environment[Island]            |14.09                    |3.65 |3.86  |< .001 |
# |Group[KLE]                     |−0.07                    |1.82 |−0.04 |0.968  |
# |Group[MLE]                     |1.95                     |1.53 |1.27  |0.204  |
# |Environment[Island]:Group[KLE] |6.61                     |4.80 |1.38  |0.168  |
# |Environment[Island]:Group[KLE] |−2.58                    |4.82 |−0.53 |0.593  |
# |Environment[Island]:Group[MLE] |7.03                     |3.86 |1.82  |0.069  |
# |Environment[Island]:Group[MLE] |−4.22                    |3.84 |−1.10 |0.271  |

# pairwise comparisons ...

pairwise <- model %>%
  emmeans(~ environment * group) %>%
  contrast('revpairwise', simple = 'environment', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# output:
# contrast       group estimate    SE  df z.ratio p.value
# long - short   ENS      4.801 4.639 Inf   1.035  0.8618
# island - short ENS     16.362 4.585 Inf   3.569  0.0023 **
# island - long  ENS     11.561 1.361 Inf   8.497  <.0001 ***
# long - short   KLE     11.407 5.171 Inf   2.206  0.1523
# island - short KLE     13.784 5.184 Inf   2.659  0.0477 *
# island - long  KLE      2.377 0.711 Inf   3.346  0.0053 **
# long - short   MLE     11.833 3.136 Inf   3.773  0.0011 **
# island - short MLE     12.138 3.162 Inf   3.839  0.0008 ***
# island - long  MLE      0.305 0.547 Inf   0.558  0.9921

# table ...

dt <- pairwise %>%
  select(group, contrast, estimate, SE, z.ratio, p.value) %>%
  my_table()

dt %>% write_rds('tables/src_ept_resumption_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ept_resumption_posthoc.csv')
dt <- read_rds('tables/src_ept_resumption_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Contrast       |Estimate~($\hat{\beta}$) |$SE$ |$z$  |$p$    |
# |:-----|:--------------|:------------------------|:----|:----|:------|
# |ENS   |Long − Short   |4.80                     |4.64 |1.03 |0.862  |
# |ENS   |Island − Short |16.36                    |4.58 |3.57 |0.002  |
# |ENS   |Island − Long  |11.56                    |1.36 |8.50 |< .001 |
# |KLE   |Long − Short   |11.41                    |5.17 |2.21 |0.152  |
# |KLE   |Island − Short |13.78                    |5.18 |2.66 |0.048  |
# |KLE   |Island − Long  |2.38                     |0.71 |3.35 |0.005  |
# |MLE   |Long − Short   |11.83                    |3.14 |3.77 |0.001  |
# |MLE   |Island − Short |12.14                    |3.16 |3.84 |< .001 |
# |MLE   |Island − Long  |0.31                     |0.55 |0.56 |0.992  |

# ---------------------------------------------------------------------------- #
# plot: token proficiency ####
# ---------------------------------------------------------------------------- #

# check participants with no proficiency scores ...

check <- ep %>%
  group_by(study, group, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  filter(is.na(proficiency) == TRUE)

# prep data for orc modeling ...

ep1 <- ep %>%
  filter(is.na(proficiency) == FALSE,
         type %in% c('gap', 'resumption', 'other'),
         environment %in% c('short', 'long', 'island')) %>%
  filter(group %in% c('KLE', 'MLE')) %>%
  group_by(study, group, participant, proficiency, environment, type) %>%
  count(type) %>%
  ungroup() %>%
  complete(type, nesting(study, group, participant, proficiency, environment), fill = list(n = 0)) %>%
  select(study, group, participant, proficiency, environment, type, n) %>%
  mutate_at(c('study', 'group', 'participant', 'environment', 'type'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'),
         type = fct_relevel(type, 'other', 'gap', 'resumption')) %>%
  arrange(study, group, participant, proficiency, environment)

temp <- ep1 %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

ep1 <- ep1 %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(environment = fct_recode(environment, 
                                  'Short' = 'short', 
                                  'Long' = 'long', 
                                  'Island' = 'island'))

summary(ep1)

p1 <- ggplot(data=filter(ep1, study == 'ORC'))
p2 <- ggplot(data=filter(ep1, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = n, col = type, fill = type),
  # geom_hline(yintercept = 2.5),
  # geom_vline(xintercept = 0),
  geom_smooth(aes(linetype = type), method=lm, alpha = .75, linewidth = .5),
  #geom_point(shape = 1, alpha = .5),
  theme_classic(),
  scale_x_continuous(name = 'Proficiency (*z*-score)', breaks = c(-2, 0, 2)),
  scale_y_continuous(name = 'Tokens', breaks = c(0, 2.5, 5)),
  coord_cartesian(ylim = c(-.1, 5.1), xlim = c(-2, 2)),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000', 'coral'),
                      breaks = c('gap', 'resumption', 'other'),
                      labels = c('Gapping', 'Resumption', 'Other')),
  scale_fill_manual(name='Dependency', 
                      values=c('#c1d2ff', '#ffdf99', '#ffbfa7'),
                      breaks = c('gap', 'resumption', 'other'),
                      labels = c('Gapping', 'Resumption', 'Other')),
  scale_linetype_manual(name = 'Dependency',
                        values = c('41', 'solid', '7222'),
                        breaks = c('gap', 'resumption', 'other'),
                        labels = c('Gapping', 'Resumption', 'Other')),
  theme(text = element_text(size = 12),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right",
        axis.title.x = element_markdown(),
        strip.text.y = element_markdown()),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ept_token_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_ept_token_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# model: resumption proficiency ####
# ---------------------------------------------------------------------------- #

# check participants with no proficiency scores ...

check <- ep %>%
  group_by(study, group, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  filter(is.na(proficiency) == TRUE)

# prep data for orc modeling ...

md1 <- ep %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(resumption = case_when(type == 'resumption' ~ TRUE,
                                TRUE ~ FALSE)) %>%
  filter(group %in% c('KLE', 'MLE')) %>%
  mutate(group = fct_drop(group)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  select(study, group, participant, item, environment, resumption, type, proficiency) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'environment', 'type'), factor) %>%
  filter(is.na(environment) == FALSE) %>%
  mutate(environment = fct_drop(environment))

temp <- md1 %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md1 <- md1 %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md1)

contrasts(md1$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md1$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md1$environment)
contrasts(md1$group)

# prep data for src modeling ...

md2 <- md1 %>%
  filter(environment %in% c('long', 'island')) %>%
  mutate(environment = fct_drop(environment))

summary(md2)

contrasts(md2$environment) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md2$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md2$environment)
contrasts(md2$group)

# plot ...

plot1 <- md1 %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant, proficiency, environment) %>%
  summarise(resumption = mean(resumption, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

plot2 <- md2 %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant, proficiency, environment) %>%
  summarise(resumption = mean(resumption, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot1, study == 'ORC'))
p2 <- ggplot(data=filter(plot2, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = resumption * 100),
  geom_hline(yintercept = 50),
  geom_vline(xintercept = 0),
  geom_smooth(method=lm, col = '#ffb000', fill = '#ffb000'),
  #geom_point(shape = 1, alpha = .5),
  theme_classic(),
  scale_x_continuous(name='Proficiency (z-score)'),
  scale_y_continuous(name="Resumption (%)", breaks = c(0, 50, 100)),
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = .5), 
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(-20, 120))
ggsave('plots/orc_ept_resumption_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(-20, 120))
ggsave('plots/src_ept_resumption_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model <- glmer(resumption ~ environment * group * proficiency + 
                  (1 + environment | participant) + 
                  (1 + environment + group | item),
                data = filter(md1, study == 'ORC'),
                family = binomial,
                control = glmerControl(optimizer = 'bobyqa', 
                                       optCtrl = list(maxfun = 1e6)))
summary(model)
toc()
beep()

# ... 431.06 sec elapsed 

# save ...

model %>% write_rds('models/orc_ept_resumption_proficiency_model.rds')
model <- read_rds('models/orc_ept_resumption_proficiency_model.rds')

# model summary ...

summary(model)

# Parameter                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                      -8.7188     1.7873  -4.878 1.07e-06 ***
# environment2                      8.8075     4.8980   1.798   0.0722 .  
# environment3                      7.8948     4.9452   1.596   0.1104    
# group2                            0.4181     1.2093   0.346   0.7295    
# proficiency                      -1.7096     0.7307  -2.340   0.0193 *  
# environment2:group2              -1.4956     2.3267  -0.643   0.5204    
# environment3:group2              -2.2347     2.5005  -0.894   0.3715    
# environment2:proficiency          2.7883     1.6575   1.682   0.0925 .  
# environment3:proficiency          2.5972     1.7237   1.507   0.1319    
# group2:proficiency                1.1542     1.2470   0.926   0.3547    
# environment2:group2:proficiency  -1.9524     2.5323  -0.771   0.4407    
# environment3:group2:proficiency  -2.2889     2.6366  -0.868   0.3853 

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3')) %>%
  my_table()

dt %>% kbl(format = 'markdown')

# |Parameter                                  |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:------------------------------------------|:------------------------|:----|:-----|:------|
# |(Intercept)                                |−8.72                    |1.79 |−4.88 |< .001 |
# |Environment[Long]                          |8.81                     |4.90 |1.80  |0.072  |
# |Environment[Island]                        |7.89                     |4.95 |1.60  |0.110  |
# |Group[MLE]                                 |0.42                     |1.21 |0.35  |0.730  |
# |Proficiency                                |−1.71                    |0.73 |−2.34 |0.019  |
# |Environment[Long]:Group[MLE]               |−1.50                    |2.33 |−0.64 |0.520  |
# |Environment[Island]:Group[MLE]             |−2.23                    |2.50 |−0.89 |0.371  |
# |Environment[Long]:Proficiency              |2.79                     |1.66 |1.68  |0.093  |
# |Environment[Island]:Proficiency            |2.60                     |1.72 |1.51  |0.132  |
# |Group[MLE]:Proficiency                     |1.15                     |1.25 |0.93  |0.355  |
# |Environment[Long]:Group[MLE]:Proficiency   |−1.95                    |2.53 |−0.77 |0.441  |
# |Environment[Island]:Group[MLE]:Proficiency |−2.29                    |2.64 |−0.87 |0.385  |

dt %>% write_rds('tables/orc_ept_resumption_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ept_resumption_proficiency_model.csv')
dt <- read_rds('tables/orc_ept_resumption_proficiency_model.rds')

# joint tests ...

joint_tests(model)

# model term                    df1 df2 F.ratio p.value
# environment                     2 Inf   1.833  0.1599
# group                           1 Inf   0.120  0.7295
# proficiency                     1 Inf   5.474  0.0193
# environment:group               2 Inf   0.522  0.5935
# environment:proficiency         2 Inf   1.479  0.2278
# group:proficiency               1 Inf   0.857  0.3547
# environment:group:proficiency   2 Inf   0.381  0.6828

# posthoc tests ...

trends <- emtrends(model, ~ group * environment, var = 'proficiency', adjust = 'mvt')
test <- test(trends)
test

# group environment proficiency.trend    SE  df z.ratio p.value
# KLE   short                  -4.789 2.422 Inf  -1.977  0.2355
# MLE   short                  -2.221 2.037 Inf  -1.090  0.8217
# KLE   long                   -1.024 0.589 Inf  -1.738  0.3702
# MLE   long                   -0.409 0.528 Inf  -0.775  0.9566
# KLE   island                 -1.047 0.817 Inf  -1.282  0.6953
# MLE   island                 -0.768 0.769 Inf  -0.999  0.8717

# table ...

dt <- test %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()
 
dt %>% kbl(format = 'markdown')

# |Group |Environment | Proficiency Trend|   SE|     z|     p|
# |:-----|:-----------|-----------------:|----:|-----:|-----:|
# |KLE   |Short       |             -4.79| 2.42| -1.98| 0.235|
# |KLE   |Long        |             -1.02| 0.59| -1.74| 0.370|
# |KLE   |Island      |             -1.05| 0.82| -1.28| 0.695|
# |MLE   |Short       |             -2.22| 2.04| -1.09| 0.822|
# |MLE   |Long        |             -0.41| 0.53| -0.77| 0.957|
# |MLE   |Island      |             -0.77| 0.77| -1.00| 0.872|

dt %>% write_rds('tables/orc_ept_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ept_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_ept_proficiency_posthoc.rds')

# plot ...

temp <- md1 %>%
  filter(study == 'ORC') %>%
  group_by(group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

plot <- emmip(model, environment ~ proficiency | group, cov.reduce = range, 
              CIs = TRUE, plotit = FALSE) %>%
  left_join(temp, by = 'group') %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')'))

ggplot(plot, aes(x = xvar, y = yvar, group = environment, col = environment, linetype = environment)) +
  theme_classic() +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = environment), alpha = .1, col = NA) +
  geom_line(lwd = .75) +
  scale_y_continuous(name='Logits') +
  scale_x_continuous(name='Proficiency (z-score)') +
  scale_linetype_manual(name = 'Environment', 
                        values = c('longdash', 'solid', 'twodash'),
                        labels = c('Short', 'Long', 'Island')) +
  scale_color_manual(name = 'Environment',
                     values = c('#619CFF', '#00BA38', '#F8766D'),
                     labels = c('Short', 'Long', 'Island')) +
  scale_fill_manual(name = 'Environment',
                    values = c('#619CFF', '#00BA38', '#F8766D'),
                    labels = c('Short', 'Long', 'Island')) +
  theme(legend.key.size = unit(2,"line")) +
  facet_wrap(~group_label)

ggsave("plots/orc_ept_resumption_proficiency_posthoc.png", width=6.5, height=2, dpi=600)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# load SRC elicited production task data ...

src <- read_csv('data/ept_modeling.csv') %>%
  filter(study == 'SRC')

# add dummy participants (with only resumption responses) to ENS and KLE groups ...

dummy <- read_csv('data/ept_dummy.csv') %>% 
  filter(group %in% c('ENS', 'KLE'))

md <- bind_rows(src, dummy) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'environment', 'type'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(study = fct_drop(study)) %>%
  mutate(environment = fct_drop(environment))

# add proficiency scores ...

prof <- read_csv("data/proficiency.csv")

prof %>% write_csv('data/proficiency.csv')

md <- md %>%
  left_join(prof)

md <- md %>%
  mutate_at(c('study', 'group', 'participant'), factor)

# inspect dataframe ...

summary(md)

# set contrasts ...

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$environment)
contrasts(md$group)

# model: proficiency effects for resumption rates

md <- md %>%
  filter(is.na(proficiency) == FALSE) %>%
  filter(group %in% c('KLE', 'MLE')) %>%
  mutate(group = fct_drop(group)) %>%
  filter(is.na(environment) == FALSE) %>%
  mutate(environment = fct_drop(environment))

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$environment)
contrasts(md$group)

# ... src model (final)

# fit model ...

tic()
model1 <- glmer(resumption ~ environment * group * proficiency +
                  (1 + environment | participant) + 
                  (1 | item), 
                data = md, family = binomial, 
                control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
summary(model1)
toc()
beep()

# 61.02 sec elapsed 
# no warnings
# nothing is significant

# Estimate Std. Error z value Pr(>|z|)
# (Intercept)                      -7.8753     9.5743  -0.823    0.411
# environment2                     19.0551    28.6927   0.664    0.507
# environment3                     19.9757    28.6939   0.696    0.486
# group2                            8.9218    19.1409   0.466    0.641
# proficiency                       0.7335     8.4707   0.087    0.931
# environment2:group2             -17.2271    57.3378  -0.300    0.764
# environment3:group2             -18.7813    57.3412  -0.328    0.743
# environment2:proficiency          0.5827    25.3853   0.023    0.982
# environment3:proficiency          1.4872    25.3859   0.059    0.953
# group2:proficiency               -0.4140    16.9321  -0.024    0.980
# environment2:group2:proficiency   2.6451    50.7100   0.052    0.958
# environment3:group2:proficiency   0.1379    50.7093   0.003    0.998

# table ...

dt <- summary(model1)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3')) %>%
  my_table()

dt %>% write_rds('tables/src_ept_resumption_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ept_resumption_proficiency_model.csv')
dt <- read_rds('tables/src_ept_resumption_proficiency_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                  |Estimate~($\hat{\beta}$) |$SE$  |$z$   |$p$   |
# |:------------------------------------------|:------------------------|:-----|:-----|:-----|
# |(Intercept)                                |−7.88                    |9.57  |−0.82 |0.411 |
# |Environment[Long]                          |19.06                    |28.69 |0.66  |0.507 |
# |Environment[Island]                        |19.98                    |28.69 |0.70  |0.486 |
# |Group[MLE]                                 |8.92                     |19.14 |0.47  |0.641 |
# |Proficiency                                |0.73                     |8.47  |0.09  |0.931 |
# |Environment[Long]:Group[MLE]               |−17.23                   |57.34 |−0.30 |0.764 |
# |Environment[Island]:Group[MLE]             |−18.78                   |57.34 |−0.33 |0.743 |
# |Environment[Long]:Proficiency              |0.58                     |25.39 |0.02  |0.982 |
# |Environment[Island]:Proficiency            |1.49                     |25.39 |0.06  |0.953 |
# |Group[MLE]:Proficiency                     |−0.41                    |16.93 |−0.02 |0.980 |
# |Environment[Long]:Group[MLE]:Proficiency   |2.65                     |50.71 |0.05  |0.958 |
# |Environment[Island]:Group[MLE]:Proficiency |0.14                     |50.71 |0.00  |0.998 |

# posthoc tests ...

trends <- emtrends(model1, ~ group * environment, var = 'proficiency', adjust = 'mvt')
test <- test(trends)
test

# group environment proficiency.trend     SE  df z.ratio p.value
# KLE   short                  0.7143 50.718 Inf   0.014  1.0000
# MLE   short                 -0.6273  1.479 Inf  -0.424  0.9977
# KLE   long                  -0.0256  0.841 Inf  -0.030  1.0000
# MLE   long                   1.2779  0.803 Inf   1.591  0.4635
# KLE   island                 2.1326  0.864 Inf   2.468  0.0723
# MLE   island                 0.9289  0.792 Inf   1.173  0.7627

# table ...

dt <- test %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ept_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ept_proficiency_posthoc.csv')
dt <- read_rds('tables/src_ept_proficiency_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Proficiency Trend |$SE$  |$z$   |$p$    |
# |:-----|:-----------|:-----------------|:-----|:-----|:------|
# |KLE   |Short       |0.71              |50.72 |0.01  |> .999 |
# |KLE   |Long        |−0.03             |0.84  |−0.03 |> .999 |
# |KLE   |Island      |2.13              |0.86  |2.47  |0.072  |
# |MLE   |Short       |−0.63             |1.48  |−0.42 |0.998  |
# |MLE   |Long        |1.28              |0.80  |1.59  |0.464  |
# |MLE   |Island      |0.93              |0.79  |1.17  |0.763  |

# ---------------------------------------------------------------------------- #
# plot: resumption subtype by group ####
# ---------------------------------------------------------------------------- #

# remove participants who only gave nontarget responses ...

temp <- ep %>%
  filter(type != 'nontarget') %>%
  group_by(study, group, participant) %>%
  summarise() %>%
  ungroup()

ep2 <- ep %>%
  filter(participant %in% temp$participant)

check <- ep2 %>%
  filter(type == 'resumption', subtype == 'who')

# summarise for plotting ...

plot <- ep2 %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  mutate(subtype = case_when(subtype == 'who' ~ 'pronoun', TRUE ~ subtype)) %>%
  filter(type == 'resumption', is.na(environment) == FALSE) %>%
  select(study, group, group_label, participant, item, subtype) %>%
  group_by(study, group, group_label) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(study, group, group_label, subtype, total) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(prop = (count/total) * 100) %>%
  complete(subtype, 
           nesting(study, group, total), 
           fill = list(count = 0, prop = 0)) %>%
  mutate(subtype = fct_relevel(subtype, 
                               'pronoun', 
                               'fullNP', 
                               #'who'
                               )) %>%
  mutate(subtype = recode_factor(subtype, 
                                 pronoun = 'RP', 
                                 fullNP = 'Full NP', 
                                 #who = 'Who'
                                 )) %>%
  select(study, group, group_label, subtype, count, total, prop) %>%
  arrange(study, group, subtype, count, total, prop)

# define data for plots ...

p1 <- ggplot(data=filter(plot, study == 'ORC' & subtype %in% c('RP', 'Full NP')))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# styling ... 

s <- list(
  aes(x = group, y = prop, group = subtype, fill = subtype, shape = subtype, label = round(prop)),
  geom_bar(stat = "identity", col = 'black', width = .3, alpha=.8),
  geom_text(size = 3, col = "black", position = position_stack(vjust = 0.5)),
  theme_classic(),
  scale_x_discrete(name="Group"),
  scale_y_continuous(name="Responses (%)", limits=c(0, 101), breaks=c(0, 20, 40, 60, 80, 100)),
  scale_fill_manual(name="Subtype",
                    values=c('#40B0A6', '#E1BE6A', '#E54C4B')),
  theme(text = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = .5), 
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "right", 
        legend.margin=margin(1, 1, 1, 1))
)

p1 + s
ggsave('plots/orc_ept_line_subtype_group.png', width=6.5, height=2.25, dpi=600)

p2 + s
ggsave('plots/src_ept_line_subtype_group.png', width=6.5, height=2.25, dpi=600)

# table ...

table <- plot %>%
  select(-total) %>%
  rename('Group' = 'group', 
         'Count' = 'count',
         'Type' = 'subtype', 
         '%' = 'prop'
  ) %>%
  mutate(Group = as.character(Group),
         Type = as.character(Type)) %>%
  mutate_if(is.character, str_replace_all, 
            c(
              'cond1' = 'Short',
              'cond2' = 'Long',
              'cond3' = 'Island',
              'pronoun' = 'RP',
              'fullNP' = 'Full NP',
              'who' = 'Who'
            )
  )

table %>%
  filter(study == 'ORC', Type %in% c('RP', 'Full NP')) %>%
  select(-study) %>%
  kbl() %>% kable_classic()

table <- table %>% my_table()

table %>% write_rds('tables/ept_subtype.rds')
csv <- table %>% csv_table() %>% write_csv('tables/ept_subtype.csv')
table <- read_rds('tables/ept_subtype.rds')

csv_orc <- csv %>%
  filter(study == "ORC", Type != "Who") %>%
  write_csv('tables/orc_ept_subtype.csv')

csv_src <- csv %>%
  filter(study == "SRC") %>%
  write_csv('tables/src_ept_subtype.csv')

# ---------------------------------------------------------------------------- #
# plot: resumption subtype by group & environment ####
# ---------------------------------------------------------------------------- #

# remove participants who only gave nontarget responses ...

temp <- ep %>%
  filter(type != 'nontarget') %>%
  group_by(study, group, participant) %>%
  summarise() %>%
  ungroup()

ep2 <- ep %>%
  filter(participant %in% temp$participant)

check <- plot %>%
  filter(group == 'ENS', environment == 'island')

# summarise for plotting ...

plot <- ep2 %>%
  filter(type == 'resumption', is.na(environment) == FALSE) %>%
  select(study, group, participant, item, environment, subtype) %>%
  group_by(study, group, environment) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(study, group, environment, subtype, total) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(prop = (count/total) * 100) %>%
  complete(subtype, 
           nesting(study, group, environment, total), 
           fill = list(count = 0, prop = 0)) %>%
  complete(environment, 
           nesting(study, group, subtype),
           fill = list(total = 0, count = 0, prop = 0)) %>%
  select(study, group, environment, subtype, count, total, prop)

# define data for plots ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# styling ...

s <- list(
  aes(x = environment, y = prop, group = subtype, col = subtype, shape = subtype),
  geom_hline(yintercept = 50),
  geom_line(position = position_dodge(width = .6)),
  geom_point(lwd = 3, position = position_dodge(width = .6), alpha = .5),
  theme_classic(),
  scale_y_continuous(name='Responses (%)', 
                     #limits = c(-.05, 5)
                     ),
  scale_x_discrete(name='Environment', 
                   limits = c('short', 'long', 'island'),
                   labels = c('Short', 'Long', 'Island')
  ),
  scale_colour_manual(name='Subtype', 
                      values=c(
                        '#648fff', 
                                 'gray80', 
                                 #'gray60', 
                                 '#ffb000'),
                                 breaks = c(
                                   'pronoun', 
                                   'fullNP', 
                                   'who'
                                 ),
                      labels = c(
                        'RP', 
                        'Full NP', 
                        'Who'
                      )
  ),
  scale_shape_manual(name='Subtype', 
                     values=c(
                       16, 
                       17, 
                       #18, 
                       15
                     ),
                     breaks = c(
                       'pronoun', 
                       'fullNP', 
                       'who'
                     ),
                     labels = c(
                       'RP', 
                       'Full NP', 
                       'Who'
                     )
  ),
  theme(text = element_text(size = 12),
        legend.position = 'right',
        legend.margin = margin(t = -.4, unit = 'cm'),
        plot.title = element_text(size = 12, hjust = .5)),
  facet_manual(~group,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ept_line_subtype.png', width=6.5, height=2.25, dpi=600)

p2 + s
ggsave('plots/src_ept_line_subtype.png', width=6.5, height=2.25, dpi=600)

# table ...

table <- plot %>%
  select(-group_label) %>%
  rename('Group' = 'group', 
         'Environment' = 'condition',
         'Subtype' = 'subtype', 
         'Mean' = 'mean',
         'LCL' = 'lcl',
         'UCL' = 'ucl',
         '%' = 'prop'
  ) %>%
  mutate(Group = as.character(Group),
         Environment = as.character(Environment),
         Subtype = as.character(Subtype)) %>%
  mutate_if(is.character, str_replace_all, 
            c(
              'cond1' = 'Short',
              'cond2' = 'Long',
              'cond3' = 'Island',
              'pronoun' = 'RP',
              'fullNP' = 'Full NP',
              'who' = 'Who'
            )
  )

table %>% write_rds('tables/ept_subtype.rds')
table %>% write_csv('tables/ept_subtype.csv')

# ============================================================================ #
# ** SELF-PACED READING TASK **  ####
# ============================================================================ #

# prep dataframe ...

spr <- df %>%
  filter(task == 'sprt') %>%
  arrange(study, group, participant) %>%
  select_if(function(x){!all(is.na(x))}) %>%
  mutate(trial = sub('.*\\.', '', internal_node_id)) %>%
  mutate(trial = as.numeric(trial) + 1) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'condition'), factor) %>%
  mutate(rt = as.numeric(as.character(rt))) %>%
  mutate(class = factor(case_when(condition %in% c('grammatical', 'ungrammatical') ~ 'filler',
                                  condition %in% c('cond1', 'cond2', 'cond3', 'cond4', 'cond5', 'cond6') ~ 'critical'))) %>%
  mutate(set = case_when(item %in% c('item31', 'item32', 'item33', 'item34', 'item35', 'item36', 'item37') ~ 1,
                         item %in% c('item38', 'item39', 'item40', 'item41', 'item42', 'item43', 'item44') ~ 1,
                         item %in% c('item45', 'item46', 'item47', 'item48', 'item49', 'item50', 'item51') ~ 2,
                         item %in% c('item52', 'item53', 'item54', 'item55', 'item56', 'item57', 'item58') ~ 2,
                         item %in% c('item59', 'item60', 'item61', 'item62', 'item63', 'item64', 'item65') ~ 3,
                         item %in% c('item66', 'item67', 'item68', 'item69', 'item70', 'item71', 'item72') ~ 4))

check1 <- df %>%
  group_by(study, group) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

check2 <- spr %>%
  group_by(study, group) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

# participant exclusions ...

temp <- spr %>%  
  group_by(study, participant, class) %>%
  summarise(acc_rate = mean(as.logical(accuracy))) %>%
  ungroup() %>%
  group_by(study, class) %>%
  mutate(mean = mean(acc_rate),
         sd = sd(acc_rate)) %>%
  ungroup() %>%
  mutate(cutoff = mean - 3 * sd) %>%
  mutate(outlier = case_when(acc_rate < cutoff ~ TRUE,
                             TRUE ~ FALSE)) %>%
  group_by(study, participant) %>%
  summarise(exclude = case_when(mean(outlier) == 1 ~ TRUE, TRUE ~ FALSE),
            acc_crit = acc_rate[class == 'critical'],
            cutoff_crit = mean(cutoff[class == 'critical']),
            acc_fill = acc_rate[class == 'filler'],
            cutoff_fill = mean(cutoff[class == 'filler'])) %>%
  ungroup() %>%
  filter(exclude == TRUE)

spr_clean <- spr %>%
  filter(!participant %in% temp$participant)

check <- spr_clean %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

# ---------------------------------------------------------------------------- #
# data trimming ####
# ---------------------------------------------------------------------------- #

# filter to critical trials ...

spr_crit <- spr_clean %>%
  filter(!condition %in% c('grammatical', 'ungrammatical'))

# adjust region numbers ...

spr_crit <- spr_crit %>%
  mutate(region = as.numeric(region)) %>%
  mutate(region2 = case_when(study == 'ORC' ~ region - 11,
                             study == 'SRC' & environment %in% c('short', 'long') ~ region - 8,
                             study == 'SRC' & environment == 'island' ~ region - 10))

# remove extreme values ...

spr_trim <- spr_crit %>%
  filter(rt <= 6000, rt >= 100)

# calculate percent excluded ...

check <- spr_crit %>%
  mutate(keep = case_when(rt <= 6000 & rt >= 100 ~ TRUE,
                             TRUE ~ FALSE))

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# remove outliers more than 2 SDs from the mean ...

spr_rt <- spr_trim %>%
  group_by(study, group, dependency, environment, region2) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  filter(rt <= upper & rt >= lower) %>%
  select(-upper, -lower)

# calculate percent excluded ...

check <- spr_trim %>%
  group_by(study, group, dependency, environment, region2) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  mutate(keep = case_when(rt <= upper & rt >= lower ~ TRUE,
                          TRUE ~ FALSE)) %>%
  select(study, group, participant, dependency, environment, region2, rt, upper, lower, keep)

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# check distribution of raw RTs ...

check <- spr_rt %>% mutate(rt = scale(rt))

hist(check$rt[check$study == 'ORC'])
qqnorm(check$rt[check$study == 'ORC'])

hist(check$rt[check$study == 'SRC'])
qqnorm(check$rt[check$study == 'SRC'])

check <- spr_rt %>% mutate(rt = log(rt))

hist(check$rt[check$study == 'ORC'])
qqnorm(spr_rt$rt[spr_rt$study == 'ORC'])

hist(check$rt[check$study == 'SRC'])
qqnorm(check$rt[check$study == 'SRC'])

# calculate RRTs by participant for each study and add to dataframe ...

spr_rrt <- spr_rt %>%
  mutate(stimulus = str_trim(tolower(stimulus))) %>%
  mutate(length = nchar(stimulus))

spr_rrt <- spr_rrt %>%
  mutate(region = as.factor(region))

spr_rrt_orc <- spr_rrt %>%
  filter(study == 'ORC')

mod_spr_rrt_orc <- lmer(log(rt) ~ length + region + trial + (1 | participant), data = spr_rrt_orc)

spr_rrt_orc <- spr_rrt_orc %>%
  mutate(rrt = residuals(mod_spr_rrt_orc))

spr_rrt_src <- spr_rrt %>%
  filter(study == 'SRC')

mod_spr_rrt_src <- lmer(log(rt) ~ length + region + trial + (1 | participant), data = spr_rrt_src)

spr_rrt_src <- spr_rrt_src %>%
  mutate(rrt = residuals(mod_spr_rrt_src))

spr_rrt <- bind_rows(spr_rrt_orc, spr_rrt_src)

# check distribution of RRTs ...

hist(spr_rrt$rrt)
qqnorm(spr_rrt$rrt)

# ---------------------------------------------------------------------------- #
# examine effect of trial number on RT ####
# ---------------------------------------------------------------------------- #

# test whether trial number influences RTs (it does, and we need to adjust for this!)

ggplot(filter(spr_rt, study == 'ORC')) +
  aes(x=trial, y=rt, group=dependency, col=dependency, shape=dependency) +
  geom_smooth(aes(linetype = dependency), method=lm, alpha = .75, linewidth = .5) +
  theme_classic() +
  facet_grid2(vars(group), vars(environment), axes = 'all', remove_labels = FALSE)

mod <- lm(rt ~ trial + dependency + environment + group, spr_rt)
summary(mod)

ggplot(filter(spr_rrt, study == 'ORC')) +
  aes(x=trial, y=rrt, group=dependency, col=dependency, shape=dependency) +
  ylim(c(-200, 200)) +
  geom_smooth(aes(linetype = dependency), method=lm, alpha = .75, linewidth = .5) +
  theme_classic() +
  facet_grid2(vars(group), vars(environment), axes = 'all', remove_labels = FALSE)

mod <- lm(rrt ~ trial + dependency + environment + group, spr_rrt)
summary(mod)

# ---------------------------------------------------------------------------- #
# plot: initial inspection ####
# ---------------------------------------------------------------------------- #

spr_rt_test <- spr_crit %>%
  filter(rt <= 6000, rt >= 100) %>%
  group_by(study, group, dependency, environment, region2) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  filter(rt <= upper & rt >= lower)
  # mutate(rt = case_when(rt > upper ~ upper,
  #                       rt < lower ~ lower,
  #                       TRUE ~ rt))

p1 <- ggplot(data = filter(spr_rt_test, study == 'ORC'))
p2 <- ggplot(data = filter(spr_rt_test, study == 'SRC'))

s <- list(
  aes(x = environment, y = rt, col = dependency),
  # geom_violin(position = position_dodge(width = .5),
  # alpha = .5,
  # col = NA
  # ),
  geom_boxplot(outlier.size = .5,
               width = .25,
               position = position_dodge(width = .5, preserve = 'single')
  ),
  # geom_jitter(),
  theme_classic(),
  scale_x_discrete(name = 'environment',
                   limits = c('short', 'long', 'island')),
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 9),
        legend.margin = margin(-10, 0, 0, 0)),
  facet_wrap(~group) 
)

p1 + s 
p2 + s 

min(spr_rt_test$rt)

ggplot(filter(spr_rt_test, study == 'ORC')) +
  aes(x=environment, y=log(rt), group=item) +
  geom_smooth(method=lm, alpha = .5, linewidth = .1, se = FALSE) +
  theme_classic() 

# ---------------------------------------------------------------------------- #
# plot: full-sentence raw RTs ####
# ---------------------------------------------------------------------------- #

temp <- spr_crit %>%
  group_by(study, group, participant, item, environment, dependency) %>%
  summarise(rt = sum(rt))

temp1 <- temp %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')'))

# remove extreme values ...

temp2 <- temp1 %>%
  filter(rt <= 30000, rt >= 100)

# calculate percent excluded ...

check <- temp1 %>%
  mutate(keep = case_when(rt <= 30000 & rt >= 100 ~ TRUE,
                          TRUE ~ FALSE))

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# remove outliers more than 2 SDs from the mean ...

temp3 <- temp2 %>%
  group_by(study, group, dependency, environment) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  filter(rt <= upper & rt >= lower) %>%
  select(-upper, -lower)

# calculate percent excluded ...

check <- temp2 %>%
  group_by(study, group, dependency, environment) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  mutate(keep = case_when(rt <= upper & rt >= lower ~ TRUE,
                          TRUE ~ FALSE)) %>%
  select(study, group, participant, dependency, environment, 
         rt, upper, lower, keep)

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# define rt dataframe

srt <- temp3

# plot ...

plot <- srt %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean = mean(rt, na.rm=T),
            lcl = my.lcl(rt),
            ucl = my.ucl(rt)) %>%
  ungroup() %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island"))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean/1000, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl/1000, ymax=ucl/1000), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Environment",
  ),
  scale_y_continuous(name="RT (s)",
                     breaks = c(4, 6, 8, 10)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~group_label, design = 'ABC', scales = 'free_x')
) 

p1 + s
ggsave("plots/orc_ajt_srt.png", width=6.5, height=2.25, dpi=600)

p2 + s
ggsave("plots/src_ajt_srt.png", width=6.5, height=2.25, dpi=600)

# table...

dt <- plot %>%
  select(study, pairing, dependency, environment, mean, lcl, ucl) %>%
  my_table() %>%
  csv_table()

dt %>% write_csv("tables/spr_srt.csv")

# ---------------------------------------------------------------------------- #
# plot: raw RTs ####
# ---------------------------------------------------------------------------- #

# summarise ...
  
plot <- spr_rt %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(environment = as.factor(environment)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, region2, dependency, environment, condition) %>%
  summarise(mean_rt = mean(rt, na.rm=T),
            lcl = my.lcl(rt),
            ucl = my.ucl(rt)) %>%
  ungroup() %>%
  filter(region2 %in% c(-1, 0, 1, 2, 3, 4)) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island')) %>%
  group_by(region2) %>%
  mutate(width = 0.01 * n()) %>%
  ungroup()

# plot ...

a1 <- ggplot(data=filter(plot, study == 'ORC'))
a2 <- ggplot(data=filter(plot, study == 'SRC'))

b <- aes(x=region2, y=mean_rt, group=dependency, col=dependency, shape=dependency, width=width*1.5)

c1 <- annotate('rect', xmin = 0.5, xmax = 3.5, 
               ymin = min(plot$lcl[plot$study == 'ORC']) - 5, 
               ymax = max(plot$ucl[plot$study == 'ORC']) + 5, 
               alpha = 0, col = 'black', linetype = 1)

c2 <- annotate('rect', xmin = 0.5, xmax = 3.5, 
               ymin = min(plot$lcl[plot$study == 'SRC']) - 10, 
               ymax = max(plot$ucl[plot$study == 'SRC']) + 10, 
               alpha = 0, col = 'black', linetype = 1)

d <- list(
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl, col = dependency), 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1,
             position = position_dodge(width = .4),
             shape = 21,
             fill = 'white',
             show.legend = FALSE),
  theme_classic(),
  scale_y_continuous(name='RT (ms)'),
  scale_x_continuous(name='Region', limits=c(-1.25, 4.25), breaks = c(-1, 0, 1, 2, 3, 4)),
  scale_colour_manual(name='Dependency:', values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        legend.position = 'bottom',
        plot.margin = margin(2, 2, 0, 2),
        strip.text.y = element_markdown(),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(vars(group_label), vars(environment), axes = 'all', remove_labels = FALSE)
)

a1 + b + c1 + d
ggsave('plots/orc_spr_rt.png', width=6.5, height=4.5, dpi=600)

a2 + b + c2 + d
ggsave('plots/src_spr_rt.png', width=6.5, height=4.5, dpi=600)

# table ...

dt <- plot %>%
  select(study, group, region2, environment, 
         dependency, mean_rt, lcl, ucl) %>%
  arrange(study, group, region2, environment, dependency) %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  my_table()

dt %>% write_rds('tables/spr_rt.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_rt.csv')
dt <- read_rds('tables/spr_rt.rds')

# ---------------------------------------------------------------------------- #
# plot: log RTs ####
# ---------------------------------------------------------------------------- #

# summarise ...

plot <- spr_rt %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(environment = as.factor(environment)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, region2, dependency, environment, condition) %>%
  summarise(mean_rt = mean(log(rt), na.rm=T),
            lcl = my.lcl(log(rt)),
            ucl = my.ucl(log(rt))) %>%
  ungroup() %>%
  filter(region2 %in% c(-1, 0, 1, 2, 3, 4)) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island')) %>%
  group_by(region2) %>%
  mutate(width = 0.01 * n()) %>%
  ungroup()

# plot ...

a1 <- ggplot(data=filter(plot, study == 'ORC'))
a2 <- ggplot(data=filter(plot, study == 'SRC'))

b <- aes(x=region2, y=mean_rt, group=dependency, col=dependency, shape=dependency, width=width*1.5)

c1 <- annotate('rect', xmin = 0.5, xmax = 3.5, 
               ymin = min(plot$lcl[plot$study == 'ORC']) - .02, 
                ymax = max(plot$ucl[plot$study == 'ORC']) + .02, 
                alpha = 0, col = 'black', linetype = 1)

c2 <- annotate('rect', xmin = 0.5, xmax = 3.5, 
               ymin = min(plot$lcl[plot$study == 'SRC']) - .02, 
               ymax = max(plot$ucl[plot$study == 'SRC']) + .02, 
               alpha = 0, col = 'black', linetype = 1)

d <- list(
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl, col = dependency), 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1,
             position = position_dodge(width = .4),
             shape = 21,
             fill = 'white',
             show.legend = FALSE),
  theme_classic(),
  scale_y_continuous(name='Log RT (ms)'),
  scale_x_continuous(name='Region', limits=c(-1.25, 4.25), breaks = c(-1, 0, 1, 2, 3, 4)),
  scale_colour_manual(name='Dependency:', values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        legend.position = 'bottom',
        plot.margin = margin(2, 2, 0, 2),
        strip.text.y = element_markdown(),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(vars(group_label), vars(environment), axes = 'all', remove_labels = FALSE)
)

a1 + b + c1 + d
ggsave('plots/orc_spr_logrt.png', width=6.5, height=4, dpi=600)

a2 + b + c2 + d
ggsave('plots/src_spr_logrt.png', width=6.5, height=4, dpi=600)

# table ...

dt <- plot %>%
  select(study, group, region2, environment, 
         dependency, mean_rt, lcl, ucl) %>%
  arrange(study, group, region2, environment, dependency) %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  my_table()

dt %>% write_rds('tables/spr_logrt.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_logrt.csv')
dt <- read_rds('tables/spr_logrt.rds')

# ---------------------------------------------------------------------------- #
# plot: log RTs (regions of interest only) ####
# ---------------------------------------------------------------------------- #

# summarise ...

plot <- spr_rt %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(environment = as.factor(environment)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, region2, dependency, environment, condition) %>%
  summarise(mean_rt = mean(log(rt), na.rm=T),
            lcl = my.lcl(log(rt)),
            ucl = my.ucl(log(rt))) %>%
  ungroup() %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island')) %>%
  group_by(region2) %>%
  mutate(width = 0.01 * n()) %>%
  ungroup() %>%
  mutate(region_label = case_when(region2 == 1 ~ 'Region 1',
                                  region2 == 2 ~ 'Region 2',
                                  region2 == 3 ~ 'Region 3'))

# plot ...

a1 <- ggplot(data=filter(plot, study == 'ORC'))
a2 <- ggplot(data=filter(plot, study == 'SRC'))

b <- aes(x=environment, y=mean_rt, group=dependency, col=dependency, shape=dependency, width=width*1.5)

d <- list(
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl, col = dependency), 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1,
             position = position_dodge(width = .4),
             shape = 21,
             fill = 'white',
             show.legend = FALSE),
  theme_classic(),
  scale_y_continuous(name='Log RT', 
                     breaks = c(5.5, 6, 6.5)),
  scale_x_discrete(name='Environment'),
  scale_colour_manual(name='Dependency:', values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        legend.position = 'bottom',
        plot.margin = margin(2, 2, 0, 2),
        strip.text.y = element_markdown(),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(vars(group_label), vars(region_label), axes = 'all', remove_labels = FALSE)
)

a1 + b + d + coord_cartesian(ylim = c(5.5, 6.6))
ggsave('plots/orc_spr_logrt2.png', width=6.5, height=5.5, dpi=600)

a2 + b + d + coord_cartesian(ylim = c(5.5, 6.6))
ggsave('plots/src_spr_logrt2.png', width=6.5, height=5.5, dpi=600)

# table ...

dt <- plot %>%
  select(study, group, region2, environment, 
         dependency, mean_rt, lcl, ucl) %>%
  arrange(study, group, region2, environment, dependency) %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  my_table()

dt %>% write_rds('tables/spr_logrt2.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_logrt2.csv')

# ---------------------------------------------------------------------------- #
# model: log RTs ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_rt %>%
  select(study, group, participant, item, trial, region2, dependency, environment, rt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'region2', 'dependency', 'environment'), factor) %>%
  mutate(rt = as.numeric(rt)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc region 1 model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'ORC', region2 == '1'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model4)
toc()
beep()

# save ...

model4 %>% write_rds('models/orc_spr_logrt_region1_model4.rds')
model4 <- read_rds('models/orc_spr_logrt_region1_model4.rds')

# 209.52 sec elapsed

# check performance ...

check_convergence(model4) # convergence: okay
check_singularity(model4, tolerance = 1e-9) # singularity: warning
check_collinearity(model4) # independence: okay
check_heteroscedasticity(model4) # homogeneity of variance: okay
check_outliers(model4) # outliers: okay
check_autocorrelation(model4) # autocorrelation: warning
MuMIn::r.squaredGLMM(model4) # R2m = 0.05333647, R2c = 0.3134198
check_model(model4) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15) 

# table ...

dt <- summary(model4)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_region1_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_region1_model.csv')
dt <- read_rds('tables/orc_spr_logrt_region1_model.rds')

# joint tests ...

joint_tests(model4)

# post-hoc tests ...

pairwise <- model4 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_region1_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_region1_posthoc.csv')
dt <- read_rds('tables/orc_spr_logrt_region1_posthoc.rds')

# plot ...

plot(emmeans(model4, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... orc region 2 model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model6 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency * environment | participant) + 
                 (1 + dependency + environment | item), 
               data = filter(md, study == 'ORC', region2 == '2'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model6)
toc() 
beep()

# ... 327.95 sec elapsed

# check performance ...

check_convergence(model6) # convergence: okay
check_singularity(model6, tolerance = 1e-9) # singularity: warning
check_collinearity(model6) # independence: okay
check_heteroscedasticity(model6) # homogeneity of variance: warning
check_outliers(model6) # outliers: okay
check_autocorrelation(model6) # autocorrelation: warning
MuMIn::r.squaredGLMM(model6) # R2m = 0.04167573, R2c = 0.3125969
check_model(model6) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15) 

# save ...

model6 %>% write_rds('models/orc_spr_logrt_region2_model6.rds')
model6 <- read_rds('models/orc_spr_logrt_region2_model6.rds')

# summary ...

summary(model6)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)  
# (Intercept)                       5.836170   0.013252 219.571478 440.400  < 2e-16 ***
# dependency2                      -0.067218   0.010677  67.922801  -6.296 2.58e-08 ***
# environment2                      0.005533   0.011433  34.346127   0.484 0.631480    
# environment3                      0.019458   0.010605 176.664244   1.835 0.068210 .  
# group2                            0.041686   0.031738 223.784873   1.313 0.190380    
# group3                            0.153131   0.031347 223.823398   4.885 1.97e-06 ***
# dependency2:environment2         -0.048856   0.020459 296.979418  -2.388 0.017568 *  
# dependency2:environment3         -0.069789   0.021307 339.592052  -3.275 0.001164 ** 
# dependency2:group2               -0.090024   0.023953 224.849735  -3.758 0.000218 ***
# dependency2:group3               -0.096465   0.023663 225.052402  -4.077 6.34e-05 ***
# environment2:group2               0.044773   0.024701 340.838123   1.813 0.070780 .  
# environment3:group2               0.026589   0.025290 350.179457   1.051 0.293805    
# environment2:group3               0.043048   0.024404 341.167881   1.764 0.078632 .  
# environment3:group3               0.046658   0.024988 350.741678   1.867 0.062701 .  
# dependency2:environment2:group2  -0.055109   0.049288 296.370559  -1.118 0.264434    
# dependency2:environment3:group2   0.087950   0.051345 339.131186   1.713 0.087644 .  
# dependency2:environment2:group3  -0.019060   0.048693 296.703923  -0.391 0.695765    
# dependency2:environment3:group3   0.081152   0.050730 339.645514   1.600 0.110598 

# table ...

dt <- summary(model6)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_region2_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_region2_model.csv')
dt <- read_rds('tables/orc_spr_logrt_region2_model.rds')

# joint tests ...

joint_tests(model6)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  39.637  <.0001 ***
# environment                    2 Inf   1.721  0.1789
# group                          2 Inf  12.322  <.0001 ***
# dependency:environment         2 Inf   5.713  0.0033 **
# dependency:group               2 Inf  10.764  <.0001 ***
# environment:group              4 Inf   1.522  0.1928
# dependency:environment:group   4 Inf   2.354  0.0515 .

# post-hoc testing ...

pairwise <- model6 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS     0.0450 0.0243 Inf   1.857  0.4416
# pronoun - gap long        ENS     0.0209 0.0243 Inf   0.861  0.9877
# pronoun - gap island      ENS    -0.0811 0.0260 Inf  -3.122  0.0160 *
# pronoun - gap short       KLE    -0.0559 0.0280 Inf  -1.995  0.3428
# pronoun - gap long        KLE    -0.1352 0.0281 Inf  -4.813  <.0001 ***
# pronoun - gap island      KLE    -0.0941 0.0301 Inf  -3.132  0.0155 *
# pronoun - gap short       MLE    -0.0721 0.0274 Inf  -2.629  0.0742 .
# pronoun - gap long        MLE    -0.1153 0.0275 Inf  -4.198  0.0002 ***
# pronoun - gap island      MLE    -0.1171 0.0294 Inf  -3.980  0.0006 ***

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_region2_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_region2_posthoc.csv')
dt <- read_rds('tables/orc_spr_logrt_region2_posthoc.rds')

# plot ...

plot(emmeans(model6, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... orc region 3 model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency * environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'ORC', region2 == '3'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc() 
beep()

# ... 1091.78 sec elapsed

# check performance ...

check_convergence(model3) # convergence: okay
check_singularity(model3, tolerance = 1e-9) # singularity: warning
check_collinearity(model3) # independence: okay
check_heteroscedasticity(model3) # homogeneity of variance: okay
check_outliers(model3) # outliers: okay
check_autocorrelation(model3) # autocorrelation: warning
MuMIn::r.squaredGLMM(model3) # R2m = 0.03900723, R2c = 0.3633028
check_model(model3) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15) 

# save ...

model3 %>% write_rds('models/orc_spr_logrt_region3_model3.rds')
model3 <- read_rds('models/orc_spr_logrt_region3_model3.rds')

# summary ...

summary(model3)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                      5.824e+00  1.884e-02  1.796e+02 309.101  < 2e-16 ***
# dependency2                     -3.542e-02  1.060e-02  4.031e+01  -3.340 0.001814 ** 
# environment2                    -1.337e-02  1.116e-02  4.908e+01  -1.198 0.236615    
# environment3                     1.321e-02  1.132e-02  5.365e+01   1.166 0.248654    
# group2                           2.986e-03  3.985e-02  2.254e+02   0.075 0.940328    
# group3                           1.717e-01  3.968e-02  2.289e+02   4.326 2.27e-05 ***
# dependency2:environment2        -1.553e-02  2.111e-02  5.753e+02  -0.736 0.462089    
# dependency2:environment3        -8.532e-03  2.079e-02  1.308e+03  -0.410 0.681541    
# dependency2:group2              -8.702e-02  2.237e-02  2.292e+02  -3.889 0.000132 ***
# dependency2:group3              -5.644e-02  2.209e-02  2.288e+02  -2.555 0.011272 *  
# environment2:group2             -7.858e-03  2.543e-02  2.923e+02  -0.309 0.757585    
# environment3:group2              1.342e-02  2.564e-02  3.302e+02   0.523 0.601040    
# environment2:group3              2.973e-02  2.519e-02  2.940e+02   1.180 0.238943    
# environment3:group3              3.649e-02  2.539e-02  3.319e+02   1.437 0.151606    
# dependency2:environment2:group2 -7.494e-02  5.078e-02  5.738e+02  -1.476 0.140555    
# dependency2:environment3:group2  1.119e-01  5.004e-02  1.305e+03   2.235 0.025557 *  
# dependency2:environment2:group3 -1.016e-02  5.005e-02  5.720e+02  -0.203 0.839208    
# dependency2:environment3:group3  1.214e-01  4.943e-02  1.304e+03   2.456 0.014193 * 

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_region3_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_region3_model.csv')
dt <- read_rds('tables/orc_spr_logrt_region3_model.rds')

# joint tests ...

joint_tests(model3)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  11.155  0.0008
# environment                    2 Inf   2.998  0.0499
# group                          2 Inf  11.477  <.0001
# dependency:environment         2 Inf   0.272  0.7616
# dependency:group               2 Inf   8.035  0.0003
# environment:group              4 Inf   0.831  0.5052
# dependency:environment:group   4 Inf   4.168  0.0022

# post-hoc tests ...

pairwise <- model3 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS     0.0369 0.0243 Inf   1.516  0.7101
# pronoun - gap long        ENS     0.0497 0.0247 Inf   2.013  0.3319
# pronoun - gap island      ENS    -0.0494 0.0246 Inf  -2.011  0.3334
# pronoun - gap short       KLE    -0.0624 0.0282 Inf  -2.215  0.2155
# pronoun - gap long        KLE    -0.1246 0.0285 Inf  -4.369  0.0001
# pronoun - gap island      KLE    -0.0369 0.0284 Inf  -1.299  0.8542
# pronoun - gap short       MLE    -0.0566 0.0275 Inf  -2.060  0.3017
# pronoun - gap long        MLE    -0.0540 0.0279 Inf  -1.934  0.3858
# pronoun - gap island      MLE    -0.0215 0.0277 Inf  -0.776  0.9942

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_region3_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_region3_posthoc.csv')
dt <- read_rds('tables/orc_spr_logrt_region3_posthoc.rds')

# plot ...

plot(emmeans(model3, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src region 1 model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model8 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'SRC', region2 == '1'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model8)
toc() 
beep()

# ... 52.94 sec elapsed

# check performance ...

check_convergence(model8) # convergence: okay
check_singularity(model8, tolerance = 1e-9) # singularity: warning
check_collinearity(model8) # independence: okay
check_heteroscedasticity(model8) # homogeneity of variance: okay
check_outliers(model8) # outliers: okay
check_autocorrelation(model8) # autocorrelation: warning
MuMIn::r.squaredGLMM(model8) # R2m = 0.03900723, R2c = 0.3633028
check_model(model8) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15) 

# save ...

model8 %>% write_rds('models/src_spr_logrt_region1_model8.rds')
model8 <- read_rds('models/src_spr_logrt_region1_model8.rds')

# summary ...

summary(model8)

# Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                      6.218e+00  2.879e-02  1.970e+02 215.990  < 2e-16 ***
# dependency2                     -3.472e-02  1.378e-02  6.220e+01  -2.520  0.01431 *  
# environment2                    -2.377e-02  1.696e-02  3.073e+01  -1.402  0.17104    
# environment3                     3.446e-02  1.713e-02  5.683e+01   2.011  0.04905 *  
# group2                           1.878e-01  6.746e-02  1.978e+02   2.783  0.00590 ** 
# group3                           4.369e-01  6.714e-02  1.998e+02   6.507 6.04e-10 ***
# dependency2:environment2        -4.363e-02  3.075e-02  5.532e+03  -1.419  0.15603    
# dependency2:environment3        -1.301e-01  3.078e-02  5.543e+03  -4.228 2.40e-05 ***
# dependency2:group2              -3.363e-02  3.125e-02  5.560e+03  -1.076  0.28192    
# dependency2:group3              -9.794e-02  3.072e-02  5.530e+03  -3.188  0.00144 ** 
# environment2:group2             -2.411e-02  3.814e-02  5.528e+03  -0.632  0.52733    
# environment3:group2             -2.841e-02  3.815e-02  5.529e+03  -0.745  0.45638    
# environment2:group3             -9.409e-03  3.764e-02  5.542e+03  -0.250  0.80261    
# environment3:group3              3.189e-02  3.764e-02  5.540e+03   0.847  0.39701    
# dependency2:environment2:group2  4.601e-02  7.652e-02  5.527e+03   0.601  0.54771    
# dependency2:environment3:group2  1.261e-01  7.658e-02  5.546e+03   1.646  0.09973 .  
# dependency2:environment2:group3 -1.356e-01  7.523e-02  5.532e+03  -1.802  0.07154 .  
# dependency2:environment3:group3 -1.874e-02  7.523e-02  5.531e+03  -0.249  0.80331   

# table ...

dt <- summary(model8)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_region1_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_region1_model.csv')
dt <- read_rds('tables/src_spr_logrt_region1_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                     | Estimate|   SE|      t|     p|
# |:---------------------------------------------|--------:|----:|------:|-----:|
# |(Intercept)                                   |     6.22| 0.03| 215.99| 0.000| ***
# |Dependency[RP]                                |    -0.03| 0.01|  -2.52| 0.014| *
# |Environment[Long]                             |    -0.02| 0.02|  -1.40| 0.171|
# |Environment[Island]                           |     0.03| 0.02|   2.01| 0.049| *
# |Group[KLE]                                    |     0.19| 0.07|   2.78| 0.006| **
# |Group[MLE]                                    |     0.44| 0.07|   6.51| 0.000| ***
# |Dependency[RP]:Environment[Long]              |    -0.04| 0.03|  -1.42| 0.156|
# |Dependency[RP]:Environment[Island]            |    -0.13| 0.03|  -4.23| 0.000| ***
# |Dependency[RP]:Group[KLE]                     |    -0.03| 0.03|  -1.08| 0.282|
# |Dependency[RP]:Group[MLE]                     |    -0.10| 0.03|  -3.19| 0.001| **
# |Environment[Long]:Group[KLE]                  |    -0.02| 0.04|  -0.63| 0.527|
# |Environment[Island]:Group[KLE]                |    -0.03| 0.04|  -0.74| 0.456|
# |Environment[Long]:Group[MLE]                  |    -0.01| 0.04|  -0.25| 0.803|
# |Environment[Island]:Group[MLE]                |     0.03| 0.04|   0.85| 0.397|
# |Dependency[RP]:Environment[Long]:Group[KLE]   |     0.05| 0.08|   0.60| 0.548|
# |Dependency[RP]:Environment[Island]:Group[KLE] |     0.13| 0.08|   1.65| 0.100|
# |Dependency[RP]:Environment[Long]:Group[MLE]   |    -0.14| 0.08|  -1.80| 0.072| .
# |Dependency[RP]:Environment[Island]:Group[MLE] |    -0.02| 0.08|  -0.25| 0.803|

# joint tests ...

joint_tests(model8)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf   6.352  0.0117
# environment                    2 Inf   5.276  0.0051
# group                          2 Inf  21.464  <.0001
# dependency:environment         2 Inf   9.250  0.0001
# dependency:group               2 Inf   5.315  0.0049
# environment:group              4 Inf   0.819  0.5130
# dependency:environment:group   4 Inf   2.235  0.0626

# post-hoc tests ...

pairwise <- model8 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS     0.0690 0.0392 Inf   1.761  0.5192
# pronoun - gap long        ENS     0.0553 0.0391 Inf   1.412  0.7867
# pronoun - gap island      ENS    -0.0969 0.0391 Inf  -2.478  0.1129
# pronoun - gap short       KLE    -0.0220 0.0381 Inf  -0.576  0.9994
# pronoun - gap long        KLE     0.0103 0.0382 Inf   0.269  1.0000
# pronoun - gap island      KLE    -0.0618 0.0385 Inf  -1.605  0.6437
# pronoun - gap short       MLE     0.0225 0.0368 Inf   0.612  0.9991
# pronoun - gap long        MLE    -0.1268 0.0369 Inf  -3.434  0.0053 **
# pronoun - gap island      MLE    -0.1621 0.0370 Inf  -4.385  0.0001 ***

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_region1_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_region1_posthoc.csv')
dt <- read_rds('tables/src_spr_logrt_region1_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|     z|     p|
# |:-----|:-----------|:--------|--------:|----:|-----:|-----:|
# |ENS   |Short       |RP - Gap |     0.07| 0.04|  1.76| 0.519|
# |ENS   |Long        |RP - Gap |     0.06| 0.04|  1.41| 0.787|
# |ENS   |Island      |RP - Gap |    -0.10| 0.04| -2.48| 0.113|
# |KLE   |Short       |RP - Gap |    -0.02| 0.04| -0.58| 0.999|
# |KLE   |Long        |RP - Gap |     0.01| 0.04|  0.27| 1.000|
# |KLE   |Island      |RP - Gap |    -0.06| 0.04| -1.61| 0.644|
# |MLE   |Short       |RP - Gap |     0.02| 0.04|  0.61| 0.999|
# |MLE   |Long        |RP - Gap |    -0.13| 0.04| -3.43| 0.005| **
# |MLE   |Island      |RP - Gap |    -0.16| 0.04| -4.38| 0.000| ***


# plot ...

plot(emmeans(model8, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src region 2 model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency * environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'SRC', region2 == '2'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc() 
beep()

# ... 838.31 sec elapsed

# check performance ...

check_convergence(model3) # convergence: okay
check_singularity(model3, tolerance = 1e-9) # singularity: warning
check_collinearity(model3) # independence: okay
check_heteroscedasticity(model3) # homogeneity of variance: warning
check_outliers(model3) # outliers: okay
check_autocorrelation(model3) # autocorrelation: warning
MuMIn::r.squaredGLMM(model3) # R2m = 0.07347688, R2c = 0.4113944
check_model(model3) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15) 

# save ...

model3 %>% write_rds('models/src_spr_logrt_region2_model3.rds')
model3 <- read_rds('models/src_spr_logrt_region2_model3.rds')

# summary ...

summary(model3)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                       6.126852   0.021872 193.182707 280.124  < 2e-16 ***
# dependency2                      -0.081278   0.011607  98.465005  -7.003 3.12e-10 ***
# environment2                      0.034798   0.015289  45.475368   2.276 0.027606 *  
# environment3                     -0.032936   0.015764 116.006197  -2.089 0.038864 *  
# group2                            0.138743   0.053835 192.918553   2.577 0.010707 *  
# group3                            0.299089   0.053272 194.164165   5.614 6.75e-08 ***
# dependency2:environment2         -0.145891   0.028698 213.060158  -5.084 8.09e-07 ***
# dependency2:environment3         -0.113884   0.028255 193.572422  -4.031 8.00e-05 ***
# dependency2:group2               -0.024834   0.027721 382.518130  -0.896 0.370890    
# dependency2:group3               -0.054201   0.027292 380.653637  -1.986 0.047753 *  
# environment2:group2               0.065228   0.035278 187.327335   1.849 0.066035 .  
# environment3:group2              -0.057764   0.038504 186.234768  -1.500 0.135248    
# environment2:group3               0.146264   0.034809 187.886828   4.202 4.09e-05 ***
# environment3:group3              -0.008374   0.037986 186.705969  -0.220 0.825763    
# dependency2:environment2:group2  -0.099422   0.071368 212.965040  -1.393 0.165047    
# dependency2:environment3:group2   0.239677   0.070239 193.466008   3.412 0.000784 ***
# dependency2:environment2:group3  -0.125659   0.070235 212.032916  -1.789 0.075021 .  
# dependency2:environment3:group3   0.173925   0.069170 192.751265   2.514 0.012740 * 

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_region2_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_region2_model.csv')
dt <- read_rds('tables/src_spr_logrt_region2_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                     | Estimate|   SE|      t|     p|
# |:---------------------------------------------|--------:|----:|------:|-----:|
# |(Intercept)                                   |     6.13| 0.02| 280.12| 0.000| ***
# |Dependency[RP]                                |    -0.08| 0.01|  -7.00| 0.000| ***
# |Environment[Long]                             |     0.03| 0.02|   2.28| 0.028| *
# |Environment[Island]                           |    -0.03| 0.02|  -2.09| 0.039| *
# |Group[KLE]                                    |     0.14| 0.05|   2.58| 0.011| *
# |Group[MLE]                                    |     0.30| 0.05|   5.61| 0.000| ***
# |Dependency[RP]:Environment[Long]              |    -0.15| 0.03|  -5.08| 0.000| ***
# |Dependency[RP]:Environment[Island]            |    -0.11| 0.03|  -4.03| 0.000| ***
# |Dependency[RP]:Group[KLE]                     |    -0.02| 0.03|  -0.90| 0.371|
# |Dependency[RP]:Group[MLE]                     |    -0.05| 0.03|  -1.99| 0.048| *
# |Environment[Long]:Group[KLE]                  |     0.07| 0.04|   1.85| 0.066| .
# |Environment[Island]:Group[KLE]                |    -0.06| 0.04|  -1.50| 0.135|
# |Environment[Long]:Group[MLE]                  |     0.15| 0.03|   4.20| 0.000| ***
# |Environment[Island]:Group[MLE]                |    -0.01| 0.04|  -0.22| 0.826|
# |Dependency[RP]:Environment[Long]:Group[KLE]   |    -0.10| 0.07|  -1.39| 0.165|
# |Dependency[RP]:Environment[Island]:Group[KLE] |     0.24| 0.07|   3.41| 0.001| **
# |Dependency[RP]:Environment[Long]:Group[MLE]   |    -0.13| 0.07|  -1.79| 0.075| .
# |Dependency[RP]:Environment[Island]:Group[MLE] |     0.17| 0.07|   2.51| 0.013| *


# joint tests ...

joint_tests(model3)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  49.035  <.0001 ***
# environment                    2 Inf   6.720  0.0012 **
# group                          2 Inf  15.827  <.0001 ***
# dependency:environment         2 Inf  14.486  <.0001 ***
# dependency:group               2 Inf   1.985  0.1373
# environment:group              4 Inf   5.506  0.0002 ***
# dependency:environment:group   4 Inf   7.102  <.0001 ***

# post-hoc tests ...

pairwise <- model3 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS     0.0526 0.0347 Inf   1.518  0.7112
# pronoun - gap long        ENS    -0.0183 0.0378 Inf  -0.483  0.9999
# pronoun - gap island      ENS    -0.1991 0.0342 Inf  -5.815  <.0001 ***
# pronoun - gap short       KLE    -0.0190 0.0337 Inf  -0.562  0.9995
# pronoun - gap long        KLE    -0.1893 0.0369 Inf  -5.135  <.0001 ***
# pronoun - gap island      KLE    -0.0311 0.0334 Inf  -0.930  0.9799
# pronoun - gap short       MLE    -0.0177 0.0327 Inf  -0.541  0.9997
# pronoun - gap long        MLE    -0.2142 0.0356 Inf  -6.011  <.0001 ***
# pronoun - gap island      MLE    -0.0955 0.0323 Inf  -2.957  0.0276 *

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_region2_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_region2_posthoc.csv')
dt <- read_rds('tables/src_spr_logrt_region2_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|     z|     p|
# |:-----|:-----------|:--------|--------:|----:|-----:|-----:|
# |ENS   |Short       |RP - Gap |     0.05| 0.03|  1.52| 0.711|
# |ENS   |Long        |RP - Gap |    -0.02| 0.04| -0.48| 1.000|
# |ENS   |Island      |RP - Gap |    -0.20| 0.03| -5.82| 0.000| ***
# |KLE   |Short       |RP - Gap |    -0.02| 0.03| -0.56| 1.000|
# |KLE   |Long        |RP - Gap |    -0.19| 0.04| -5.13| 0.000| ***
# |KLE   |Island      |RP - Gap |    -0.03| 0.03| -0.93| 0.980|
# |MLE   |Short       |RP - Gap |    -0.02| 0.03| -0.54| 1.000|
# |MLE   |Long        |RP - Gap |    -0.21| 0.04| -6.01| 0.000| ***
# |MLE   |Island      |RP - Gap |    -0.10| 0.03| -2.96| 0.028| *


# plot ...

plot(emmeans(model3, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src region 3 model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency * environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'SRC', region2 == '3'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc() 
beep()

# 805.24 sec elapsed

# check performance ...

check_convergence(model3) # convergence: okay
check_singularity(model3, tolerance = 1e-9) # singularity: warning
check_collinearity(model3) # independence: okay
check_heteroscedasticity(model3) # homogeneity of variance: warning
check_outliers(model3) # outliers: okay
check_autocorrelation(model3) # autocorrelation: warning
MuMIn::r.squaredGLMM(model3) # R2m = 0.1163956, R2c = 0.4967666
check_model(model3) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15)

# save ...

model3 %>% write_rds('models/src_spr_logrt_region3_model3.rds') 
model3 <- read_rds('models/src_spr_logrt_region3_model3.rds')

# summary ...

summary(model3)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                       6.11657    0.02587 201.43878 236.450  < 2e-16 ***
# dependency2                      -0.03629    0.01311  70.15026  -2.767  0.00722 ** 
# environment2                      0.01885    0.01807  44.12824   1.043  0.30264    
# environment3                     -0.24375    0.02618  69.13330  -9.311 7.98e-14 ***
# group2                            0.12153    0.06121 194.57739   1.985  0.04850 *  
# group3                            0.35056    0.06128 199.06808   5.721 3.85e-08 ***
# dependency2:environment2         -0.10778    0.03338 192.97998  -3.229  0.00146 ** 
# dependency2:environment3         -0.16914    0.03099 218.19003  -5.457 1.31e-07 ***
# dependency2:group2               -0.01399    0.03025 291.82608  -0.462  0.64410    
# dependency2:group3               -0.02374    0.02977 291.38271  -0.797  0.42589    
# environment2:group2               0.02405    0.03861 193.08249   0.623  0.53402    
# environment3:group2              -0.13739    0.05358 190.29943  -2.564  0.01111 *  
# environment2:group3               0.01754    0.03814 194.02818   0.460  0.64609    
# environment3:group3              -0.26157    0.05286 190.79667  -4.948 1.64e-06 ***
# dependency2:environment2:group2  -0.15539    0.08292 192.92597  -1.874  0.06243 .  
# dependency2:environment3:group2  -0.03610    0.07699 217.68049  -0.469  0.63963    
# dependency2:environment2:group3  -0.03029    0.08156 192.28042  -0.371  0.71075    
# dependency2:environment3:group3   0.13839    0.07579 217.03398   1.826  0.06923 .  

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_region3_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_region3_model.csv')
dt <- read_rds('tables/src_spr_logrt_region3_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                     | Estimate|   SE|      t|     p|
# |:---------------------------------------------|--------:|----:|------:|-----:|
# |(Intercept)                                   |     6.12| 0.03| 236.45| 0.000| ***
# |Dependency[RP]                                |    -0.04| 0.01|  -2.77| 0.007| **
# |Environment[Long]                             |     0.02| 0.02|   1.04| 0.303|
# |Environment[Island]                           |    -0.24| 0.03|  -9.31| 0.000| ***
# |Group[KLE]                                    |     0.12| 0.06|   1.99| 0.049| *
# |Group[MLE]                                    |     0.35| 0.06|   5.72| 0.000| ***
# |Dependency[RP]:Environment[Long]              |    -0.11| 0.03|  -3.23| 0.001| **
# |Dependency[RP]:Environment[Island]            |    -0.17| 0.03|  -5.46| 0.000| ***
# |Dependency[RP]:Group[KLE]                     |    -0.01| 0.03|  -0.46| 0.644|
# |Dependency[RP]:Group[MLE]                     |    -0.02| 0.03|  -0.80| 0.426|
# |Environment[Long]:Group[KLE]                  |     0.02| 0.04|   0.62| 0.534|
# |Environment[Island]:Group[KLE]                |    -0.14| 0.05|  -2.56| 0.011| *
# |Environment[Long]:Group[MLE]                  |     0.02| 0.04|   0.46| 0.646|
# |Environment[Island]:Group[MLE]                |    -0.26| 0.05|  -4.95| 0.000| ***
# |Dependency[RP]:Environment[Long]:Group[KLE]   |    -0.16| 0.08|  -1.87| 0.062| .
# |Dependency[RP]:Environment[Island]:Group[KLE] |    -0.04| 0.08|  -0.47| 0.640|
# |Dependency[RP]:Environment[Long]:Group[MLE]   |    -0.03| 0.08|  -0.37| 0.711|
# |Dependency[RP]:Environment[Island]:Group[MLE] |     0.14| 0.08|   1.83| 0.069| .


# joint tests ...

joint_tests(model3)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf   7.658  0.0057 **
# environment                    2 Inf  62.275  <.0001 ***
# group                          2 Inf  17.104  <.0001 ***
# dependency:environment         2 Inf  14.935  <.0001 ***
# dependency:group               2 Inf   0.319  0.7266
# environment:group              4 Inf   7.764  <.0001 ***
# dependency:environment:group   4 Inf   2.586  0.0350 *

# post-hoc tests ...

pairwise <- model3 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group  estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS    0.059326 0.0410 Inf   1.447  0.7629
# pronoun - gap long        ENS    0.013437 0.0416 Inf   0.323  1.0000
# pronoun - gap island      ENS   -0.143908 0.0355 Inf  -4.052  0.0005 ***
# pronoun - gap short       KLE    0.109169 0.0401 Inf   2.721  0.0571 .
# pronoun - gap long        KLE   -0.092115 0.0404 Inf  -2.279  0.1864
# pronoun - gap island      KLE   -0.130161 0.0344 Inf  -3.779  0.0014 **
# pronoun - gap short       MLE   -0.000447 0.0386 Inf  -0.012  1.0000
# pronoun - gap long        MLE   -0.076629 0.0392 Inf  -1.955  0.3731
# pronoun - gap island      MLE   -0.065293 0.0335 Inf  -1.951  0.3760

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_region3_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_region3_posthoc.csv')
dt <- read_rds('tables/src_spr_logrt_region3_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|     z|     p|
# |:-----|:-----------|:--------|--------:|----:|-----:|-----:|
# |ENS   |Short       |RP - Gap |     0.06| 0.04|  1.45| 0.763|
# |ENS   |Long        |RP - Gap |     0.01| 0.04|  0.32| 1.000|
# |ENS   |Island      |RP - Gap |    -0.14| 0.04| -4.05| 0.000| ***
# |KLE   |Short       |RP - Gap |     0.11| 0.04|  2.72| 0.057| .
# |KLE   |Long        |RP - Gap |    -0.09| 0.04| -2.28| 0.186|
# |KLE   |Island      |RP - Gap |    -0.13| 0.03| -3.78| 0.001| **
# |MLE   |Short       |RP - Gap |     0.00| 0.04| -0.01| 1.000|
# |MLE   |Long        |RP - Gap |    -0.08| 0.04| -1.95| 0.373|
# |MLE   |Island      |RP - Gap |    -0.07| 0.03| -1.95| 0.376|


# plot ...

plot(emmeans(model3, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: log RT proficiency (single region of interest) ####
# ---------------------------------------------------------------------------- #

# check region with largest mean [RP - gap] contrast across all participants

check <- spr_rt %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  mutate(logrt = log(rt)) %>%
  group_by(study, region2) %>%
  summarise(diff = mean(logrt[dependency == 'pronoun']) - mean(logrt[dependency == 'gap'])) %>%
  ungroup()

# prep dataframe ...

md <- spr_rt %>%
  filter(group %in% c('KLE', 'MLE'), region2 %in% c(1, 2, 3)) %>%
  filter(is.na(proficiency) == FALSE) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate_at(c('environment', 'dependency', 'group'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant')) %>%
  filter(is.na(proficiency) == FALSE)

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group, group_label, participant, proficiency, dependency, environment, region2) %>%
  summarise(logrt = mean(log(rt), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC', region2 == 1))
p2 <- ggplot(data=filter(plot, study == 'SRC', region2 == 2))

s <- list(
  aes(x = proficiency, y = logrt, col = dependency, fill = dependency),
  #geom_hline(yintercept = 0),
  #geom_vline(xintercept = 0),
  geom_smooth(aes(linetype = dependency), method=lm, alpha = .75, linewidth = .5),
  #geom_point(shape = 1, alpha = .5),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)',
                     breaks = c(-2, 0, 2)),
  scale_y_continuous(name='Log RT', 
                     breaks = c(5.5, 6, 6.5)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        axis.title.x = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right",
        strip.text.y = element_markdown()),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(5.5, 6.6), xlim = c(-2, 2)) 
ggsave('plots/orc_spr_rt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(5.5, 6.6))
ggsave('plots/src_spr_rt_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(log(rt) ~ dependency * environment * group * proficiency + 
                 (1 + dependency * environment | participant) +
                 (1 + dependency + environment + group | item),
               data = filter(md, study == 'ORC', region2 == 1),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

# ... 215.78 sec elapsed
# ... singular fit

check_convergence(model3) # okay
check_singularity(model3) # warning (random effects table looks okay)
check_singularity(model3, tolerance = 1e-9) # okay
check_collinearity(model3) # independence: okay
check_heteroscedasticity(model3) # homogeneity of variance: okay
check_outliers(model3) # outliers: okay
check_autocorrelation(model3) # autocorrelation: warning
MuMIn::r.squaredGLMM(model3) # R2m = 0.04667598, R2c = 0.3131932
check_model(model3) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15)

# save ...

model3 %>% write_rds('models/orc_spr_logrt_proficiency_model.rds')
model3 <- read_rds('models/orc_spr_logrt_proficiency_model.rds')

# summary ...

summary(model3)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                   6.086006   0.018156 133.839091 335.201  < 2e-16 ***
# dependency2                                  -0.070431   0.016858  46.841052  -4.178 0.000127 ***
# environment2                                  0.012974   0.019325  34.768138   0.671 0.506425    
# environment3                                  0.019046   0.016254  60.121140   1.172 0.245919    
# group2                                        0.107705   0.036125 134.163837   2.981 0.003408 ** 
# proficiency                                  -0.043606   0.017668 131.856766  -2.468 0.014863 *  
# dependency2:environment2                     -0.060051   0.029324 186.057991  -2.048 0.041981 *  
# dependency2:environment3                      0.053315   0.028541 226.987219   1.868 0.063046 .  
# dependency2:group2                           -0.006087   0.029089 132.861736  -0.209 0.834563    
# environment2:group2                           0.011740   0.030837 132.403992   0.381 0.704041    
# environment3:group2                           0.029624   0.029580 210.288435   1.001 0.317745    
# dependency2:proficiency                      -0.012977   0.014666 134.759509  -0.885 0.377811    
# environment2:proficiency                      0.015836   0.015557 134.884722   1.018 0.310518    
# environment3:proficiency                     -0.017168   0.014886 213.288770  -1.153 0.250054    
# group2:proficiency                            0.038352   0.035332 131.799066   1.085 0.279685    
# dependency2:environment2:group2              -0.032818   0.058730 186.543939  -0.559 0.576981    
# dependency2:environment3:group2               0.038577   0.057127 227.522931   0.675 0.500177    
# dependency2:environment2:proficiency         -0.062431   0.029716 192.292055  -2.101 0.036952 *  
# dependency2:environment3:proficiency         -0.031595   0.028959 233.767538  -1.091 0.276384    
# dependency2:group2:proficiency               -0.056193   0.029360 135.082610  -1.914 0.057741 .  
# environment2:group2:proficiency              -0.080464   0.031015 133.538806  -2.594 0.010535 *  
# environment3:group2:proficiency               0.009694   0.029755 212.918472   0.326 0.744892    
# dependency2:environment2:group2:proficiency   0.066406   0.059128 188.817495   1.123 0.262820    
# dependency2:environment3:group2:proficiency  -0.030945   0.057697 231.163357  -0.536 0.592243  

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_spr_logrt_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_proficiency_model.csv')
dt <- read_rds('tables/orc_spr_logrt_proficiency_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                            |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                                          |6.09                     |0.02 |335.20 |< .001 |
# |Dependency[RP]                                       |−0.07                    |0.02 |−4.18  |< .001 |
# |Environment[Long]                                    |0.01                     |0.02 |0.67   |0.506  |
# |Environment[Island]                                  |0.02                     |0.02 |1.17   |0.246  |
# |Group[KLE]                                           |0.11                     |0.04 |2.98   |0.003  |
# |Proficiency                                          |−0.04                    |0.02 |−2.47  |0.015  |
# |Dependency[RP]:Environment[Long]                     |−0.06                    |0.03 |−2.05  |0.042  |
# |Dependency[RP]:Environment[Island]                   |0.05                     |0.03 |1.87   |0.063  |
# |Dependency[RP]:Group[KLE]                            |−0.01                    |0.03 |−0.21  |0.835  |
# |Environment[Long]:Group[KLE]                         |0.01                     |0.03 |0.38   |0.704  |
# |Environment[Island]:Group[KLE]                       |0.03                     |0.03 |1.00   |0.318  |
# |Dependency[RP]:Proficiency                           |−0.01                    |0.01 |−0.88  |0.378  |
# |Environment[Long]:Proficiency                        |0.02                     |0.02 |1.02   |0.311  |
# |Environment[Island]:Proficiency                      |−0.02                    |0.01 |−1.15  |0.250  |
# |Group[KLE]:Proficiency                               |0.04                     |0.04 |1.09   |0.280  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]               |−0.03                    |0.06 |−0.56  |0.577  |
# |Dependency[RP]:Envir.[Island]:Group[KLE]             |0.04                     |0.06 |0.68   |0.500  |
# |Dependency[RP]:Envir.[Long]:Proficiency              |−0.06                    |0.03 |−2.10  |0.037  |
# |Dependency[RP]:Envir.[Island]:Proficiency            |−0.03                    |0.03 |−1.09  |0.276  |
# |Dependency[RP]:Group[KLE]:Proficiency                |−0.06                    |0.03 |−1.91  |0.058  |
# |Environment[Long]:Group[KLE]:Proficiency             |−0.08                    |0.03 |−2.59  |0.011  |
# |Environment[Island]:Group[KLE]:Proficiency           |0.01                     |0.03 |0.33   |0.745  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]:Proficiency   |0.07                     |0.06 |1.12   |0.263  |
# |Dependency[RP]:Envir.[Island]:Group[KLE]:Proficiency |−0.03                    |0.06 |−0.54  |0.592  |

# joint tests ...

joint_tests(model3)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  17.449  <.0001 ***
# environment                                2 Inf   0.688  0.5027
# group                                      1 Inf   8.882  0.0029 **
# proficiency                                1 Inf   6.092  0.0136 *
# dependency:environment                     2 Inf   7.555  0.0005 ***
# dependency:group                           1 Inf   0.043  0.8359
# dependency:proficiency                     1 Inf   0.783  0.3762
# environment:group                          2 Inf   0.518  0.5958
# environment:proficiency                    2 Inf   2.561  0.0772 .
# group:proficiency                          1 Inf   1.178  0.2777
# dependency:environment:group               2 Inf   0.755  0.4701
# dependency:environment:proficiency         2 Inf   2.208  0.1100
# dependency:group:proficiency               1 Inf   3.663  0.0556 .
# environment:group:proficiency              2 Inf   5.442  0.0043 **
# dependency:environment:group:proficiency   2 Inf   1.429  0.2397

# post-hoc tests ...

pairwise <- model3 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       KLE     0.0524 0.0323 Inf   1.619  0.4827
# pronoun - gap long        KLE    -0.0433 0.0346 Inf  -1.251  0.7535
# pronoun - gap island      KLE     0.0362 0.0300 Inf   1.209  0.7811
# pronoun - gap short       MLE    -0.0156 0.0313 Inf  -0.500  0.9967
# pronoun - gap long        MLE    -0.0449 0.0329 Inf  -1.363  0.6742
# pronoun - gap island      MLE    -0.0627 0.0284 Inf  -2.204  0.1526

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_proficiency_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$   |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:-----|
# |KLE   |Short       |RP − Gap |0.05                     |0.03 |1.62  |0.483 |
# |KLE   |Long        |RP − Gap |−0.04                    |0.03 |−1.25 |0.754 |
# |KLE   |Island      |RP − Gap |0.04                     |0.03 |1.21  |0.781 |
# |MLE   |Short       |RP − Gap |−0.02                    |0.03 |−0.50 |0.997 |
# |MLE   |Long        |RP − Gap |−0.04                    |0.03 |−1.36 |0.674 |
# |MLE   |Island      |RP − Gap |−0.06                    |0.03 |−2.20 |0.153 |

# plot ...

emmip(model3, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4 <- lmer(log(rt) ~ dependency * environment * group * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 + dependency + environment + group | item),
               data = filter(md, study == 'SRC', region2 == 2),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model4)
toc()
beep()

# ... 75.33 sec elapsed
# ... singular fit

check_convergence(model4) # okay
check_singularity(model4, tolerance = 1e-9) # warning (random effects table looks okay)
check_collinearity(model4) # independence: okay
check_heteroscedasticity(model4) # homogeneity of variance: okay
check_outliers(model4) # outliers: okay
check_autocorrelation(model4) # autocorrelation: warning
MuMIn::r.squaredGLMM(model4) # R2m = 0.05871691, R2c = 0.315456
check_model(model4) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15)

# save ...

model4 %>% write_rds('models/src_spr_logrt_proficiency_model.rds')
model4 <- read_rds('models/src_spr_logrt_proficiency_model.rds')

# summary ...

summary(model4)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                  6.153e+00  2.214e-02  1.288e+02 277.864  < 2e-16 ***
# dependency2                                 -9.118e-02  1.452e-02  7.102e+01  -6.278 2.40e-08 ***
# environment2                                 6.147e-02  2.034e-02  4.115e+01   3.022 0.004310 ** 
# environment3                                -4.574e-02  1.814e-02  7.612e+01  -2.521 0.013783 *  
# group2                                       1.661e-01  4.385e-02  1.258e+02   3.789 0.000233 ***
# proficiency                                  2.874e-02  2.200e-02  1.285e+02   1.306 0.193775    
# dependency2:environment2                    -1.767e-01  3.276e-02  3.339e+03  -5.392 7.46e-08 ***
# dependency2:environment3                    -4.705e-02  3.267e-02  3.343e+03  -1.440 0.149863    
# dependency2:group2                          -2.688e-02  2.787e-02  1.307e+02  -0.964 0.336615    
# environment2:group2                          7.361e-02  3.580e-02  1.615e+02   2.056 0.041412 *  
# environment3:group2                          3.699e-02  3.468e-02  2.717e+02   1.067 0.287069    
# dependency2:proficiency                     -9.694e-03  1.415e-02  1.328e+02  -0.685 0.494373    
# environment2:proficiency                    -2.384e-03  1.831e-02  1.679e+02  -0.130 0.896586    
# environment3:proficiency                    -1.360e-02  1.752e-02  2.731e+02  -0.776 0.438497    
# group2:proficiency                           5.800e-02  4.401e-02  1.286e+02   1.318 0.189865    
# dependency2:environment2:group2             -3.908e-02  6.555e-02  3.359e+03  -0.596 0.551043    
# dependency2:environment3:group2             -8.303e-02  6.532e-02  3.356e+03  -1.271 0.203748    
# dependency2:environment2:proficiency        -7.251e-02  3.353e-02  3.358e+03  -2.162 0.030659 *  
# dependency2:environment3:proficiency        -5.588e-02  3.293e-02  3.365e+03  -1.697 0.089794 .  
# dependency2:group2:proficiency              -7.469e-03  2.830e-02  1.328e+02  -0.264 0.792236    
# environment2:group2:proficiency              1.413e-02  3.666e-02  1.684e+02   0.385 0.700468    
# environment3:group2:proficiency             -6.899e-03  3.509e-02  2.742e+02  -0.197 0.844293    
# dependency2:environment2:group2:proficiency  3.957e-02  6.710e-02  3.368e+03   0.590 0.555420    
# dependency2:environment3:group2:proficiency  4.076e-02  6.589e-02  3.367e+03   0.619 0.536234 

# table ...

dt <- summary(model4)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3')) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_proficiency_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_proficiency_model.csv')
dt <- read_rds('tables/src_spr_logrt_proficiency_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                            |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                                          |6.15                     |0.02 |277.86 |< .001 | ***
# |Dependency[RP]                                       |−0.09                    |0.01 |−6.28  |< .001 | ***
# |Environment[Long]                                    |0.06                     |0.02 |3.02   |0.004  | **
# |Environment[Island]                                  |−0.05                    |0.02 |−2.52  |0.014  | *
# |Group[MLE]                                           |0.17                     |0.04 |3.79   |< .001 | ***
# |Proficiency                                          |0.03                     |0.02 |1.31   |0.194  |
# |Dependency[RP]:Environment[Long]                     |−0.18                    |0.03 |−5.39  |< .001 | ***
# |Dependency[RP]:Environment[Island]                   |−0.05                    |0.03 |−1.44  |0.150  |
# |Dependency[RP]:Group[MLE]                            |−0.03                    |0.03 |−0.96  |0.337  |
# |Environment[Long]:Group[MLE]                         |0.07                     |0.04 |2.06   |0.041  | *
# |Environment[Island]:Group[MLE]                       |0.04                     |0.03 |1.07   |0.287  |
# |Dependency[RP]:Proficiency                           |−0.01                    |0.01 |−0.69  |0.494  |
# |Environment[Long]:Proficiency                        |0.00                     |0.02 |−0.13  |0.897  |
# |Environment[Island]:Proficiency                      |−0.01                    |0.02 |−0.78  |0.438  |
# |Group[MLE]:Proficiency                               |0.06                     |0.04 |1.32   |0.190  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]               |−0.04                    |0.07 |−0.60  |0.551  |
# |Dependency[RP]:Envir.[Island]:Group[MLE]             |−0.08                    |0.07 |−1.27  |0.204  |
# |Dependency[RP]:Envir.[Long]:Proficiency              |−0.07                    |0.03 |−2.16  |0.031  | *
# |Dependency[RP]:Envir.[Island]:Proficiency            |−0.06                    |0.03 |−1.70  |0.090  | .
# |Dependency[RP]:Group[MLE]:Proficiency                |−0.01                    |0.03 |−0.26  |0.792  |
# |Environment[Long]:Group[MLE]:Proficiency             |0.01                     |0.04 |0.39   |0.700  |
# |Environment[Island]:Group[MLE]:Proficiency           |−0.01                    |0.04 |−0.20  |0.844  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]:Proficiency   |0.04                     |0.07 |0.59   |0.555  |
# |Dependency[RP]:Envir.[Island]:Group[MLE]:Proficiency |0.04                     |0.07 |0.62   |0.536  |

# joint tests ...

joint_tests(model4)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  39.356  <.0001 ***
# environment                                2 Inf  12.352  <.0001 ***
# group                                      1 Inf  14.289  0.0002 ***
# proficiency                                1 Inf   1.706  0.1914
# dependency:environment                     2 Inf  15.501  <.0001 ***
# dependency:group                           1 Inf   0.927  0.3357
# dependency:proficiency                     1 Inf   0.470  0.4932
# environment:group                          2 Inf   2.198  0.1110
# environment:proficiency                    2 Inf   0.310  0.7338
# group:proficiency                          1 Inf   1.737  0.1875
# dependency:environment:group               2 Inf   0.814  0.4432
# dependency:environment:proficiency         2 Inf   2.588  0.0752 .
# dependency:group:proficiency               1 Inf   0.070  0.7918
# environment:group:proficiency              2 Inf   0.132  0.8759
# dependency:environment:group:proficiency   2 Inf   0.243  0.7839

# post-hoc tests ...

pairwise <- model4 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       KLE     0.0502 0.0351 Inf   1.432  0.6283
# pronoun - gap long        KLE    -0.0421 0.0360 Inf  -1.170  0.8102
# pronoun - gap island      KLE    -0.0260 0.0347 Inf  -0.751  0.9731
# pronoun - gap short       MLE     0.0160 0.0324 Inf   0.493  0.9971
# pronoun - gap long        MLE    -0.0367 0.0325 Inf  -1.131  0.8327
# pronoun - gap island      MLE    -0.0195 0.0314 Inf  -0.622  0.9898

# table ...

dt <- pairwise %>% #fred
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

# table <- pairwise %>%
#   rename('Group' = 'group', 'Environment' = 'environment', 
#          'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
#          'p' = 'p.value') %>%
#   select(Group, Environment, Contrast, everything()) %>%
#   mutate_at(c('Environment', 'Contrast'), as.character) %>%
#   mutate_if(is.character, str_replace_all, 
#             c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
#               'gap' = 'Gap', 'pronoun' = 'RP')) %>%
#   select(-df) %>%
#   mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
#   arrange(Group, Environment) %>%
#   mutate(across(4:6, round, 2), across(7, round, 3)) %>%
#   mutate_all(as.character())

dt %>% write_rds('tables/src_spr_logrt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_proficiency_posthoc.csv')
dt <- read_rds('tables/src_spr_logrt_proficiency_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|     z|     p|
# |:-----|:-----------|:--------|--------:|----:|-----:|-----:|
# |KLE   |Short       |RP - Gap |     0.05| 0.04|  1.43| 0.628|
# |KLE   |Long        |RP - Gap |    -0.04| 0.04| -1.17| 0.810|
# |KLE   |Island      |RP - Gap |    -0.03| 0.03| -0.75| 0.973|
# |MLE   |Short       |RP - Gap |     0.02| 0.03|  0.49| 0.997|
# |MLE   |Long        |RP - Gap |    -0.04| 0.03| -1.13| 0.833|
# |MLE   |Island      |RP - Gap |    -0.02| 0.03| -0.62| 0.990|

# plot ...

emmip(model4, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: log RT proficiency (combined region of interest) ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_rt %>%
  filter(group %in% c('KLE', 'MLE'), region2 %in% c(1, 2, 3)) %>%
  filter(is.na(proficiency) == FALSE) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate_at(c('environment', 'dependency', 'group'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, participant, proficiency, item, dependency, environment, condition) %>%
  summarise(rt = sum(rt, na.rm=T)) %>%
  ungroup() %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant')) %>%
  filter(is.na(proficiency) == FALSE)

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group, group_label, participant, proficiency, dependency, environment) %>%
  summarise(rt = mean(rt, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = rt, col = dependency, fill = dependency),
  #geom_hline(yintercept = 0),
  #geom_vline(xintercept = 0),
  geom_smooth(aes(linetype = dependency), method=lm, alpha = .75, linewidth = .5),
  #geom_point(shape = 1, alpha = .5),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name="RT (ms)"),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        axis.title.x = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right",
        strip.text.y = element_markdown()),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s #+ coord_cartesian(ylim = c(-750, 750))
ggsave('plots/orc_spr_rt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s #+ coord_cartesian(ylim = c(-500, 500))
ggsave('plots/src_spr_rt_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(log(rt) ~ dependency * environment * group * proficiency + 
                 (1 + dependency * environment | participant) +
                 (1 + dependency + environment + group | item),
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

# ... 261.41 sec elapsed

performance::check_convergence(model3) # convergence: okay
check_singularity(model3, tolerance = 1e-9) # singularity: warning
check_collinearity(model3) # independence: okay
check_heteroscedasticity(model3) # homogeneity of variance: okay
check_outliers(model3) # outliers: okay
check_autocorrelation(model3) # autocorrelation: warning
MuMIn::r.squaredGLMM(model3) # R2m = 0.05814087, R2c = 0.307952
check_model(model3) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15)
anova(model3, model4, refit = FALSE)

# save ...

model3 %>% write_rds('models/orc_spr_logrt_proficiency_model3.rds')
model3 <- read_rds('models/orc_spr_logrt_proficiency_model3.rds')

# summary ...

summary(model3)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                   6.996253   0.014991 133.467902 466.688  < 2e-16 ***
# dependency2                                  -0.068253   0.011761  32.528083  -5.803 1.81e-06 ***
# environment2                                  0.014315   0.014738  35.471778   0.971 0.337979    
# environment3                                  0.029102   0.013120  60.972811   2.218 0.030280 *  
# group2                                        0.119127   0.029248 133.387122   4.073 7.92e-05 ***
# proficiency                                  -0.028486   0.014612 132.742212  -1.950 0.053338 .  
# dependency2:environment2                     -0.081573   0.023849 140.760504  -3.420 0.000818 ***
# dependency2:environment3                      0.023948   0.022633 411.677355   1.058 0.290642    
# dependency2:group2                            0.003102   0.019852 185.046807   0.156 0.876004    
# environment2:group2                           0.007245   0.024864 167.426823   0.291 0.771115    
# environment3:group2                           0.043549   0.024640 131.763258   1.767 0.079477 .  
# dependency2:proficiency                      -0.005838   0.009994 186.209406  -0.584 0.559829    
# environment2:proficiency                      0.009877   0.012469 166.666034   0.792 0.429425    
# environment3:proficiency                     -0.010227   0.012357 131.004965  -0.828 0.409360    
# group2:proficiency                            0.025605   0.029221 132.695231   0.876 0.382483    
# dependency2:environment2:group2               0.023486   0.047709 141.039914   0.492 0.623291    
# dependency2:environment3:group2               0.025535   0.045302 412.693614   0.564 0.573283    
# dependency2:environment2:proficiency         -0.034127   0.024023 142.018409  -1.421 0.157626    
# dependency2:environment3:proficiency         -0.019349   0.022825 414.429384  -0.848 0.397094    
# dependency2:group2:proficiency               -0.052070   0.019982 185.807810  -2.606 0.009907 ** 
# environment2:group2:proficiency              -0.056783   0.024900 165.906828  -2.280 0.023855 *  
# environment3:group2:proficiency              -0.004842   0.024709 130.906434  -0.196 0.844929    
# dependency2:environment2:group2:proficiency   0.070351   0.047873 140.350867   1.470 0.143926    
# dependency2:environment3:group2:proficiency   0.010366   0.045514 411.959996   0.228 0.819958 

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_spr_logrt_proficiency_model.rds')
dt <- read_rds('tables/orc_spr_logrt_proficiency_model.rds')

class(table$Parameter)

# joint tests ...

joint_tests(model3)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  33.684  <.0001 ***
# environment                                2 Inf   2.548  0.0782 .
# group                                      1 Inf  16.594  <.0001 ***
# proficiency                                1 Inf   3.801  0.0512 .
# dependency:environment                     2 Inf  10.324  <.0001 ***
# dependency:group                           1 Inf   0.024  0.8773
# dependency:proficiency                     1 Inf   0.341  0.5591
# environment:group                          2 Inf   1.820  0.1619
# environment:proficiency                    2 Inf   1.347  0.2601
# group:proficiency                          1 Inf   0.768  0.3809
# dependency:environment:group               2 Inf   0.193  0.8244
# dependency:environment:proficiency         2 Inf   1.030  0.3570
# dependency:group:proficiency               1 Inf   6.791  0.0092 **
# environment:group:proficiency              2 Inf   3.239  0.0392 *
# dependency:environment:group:proficiency   2 Inf   1.213  0.2972

# post-hoc tests ...

pairwise <- model3 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       KLE     0.0515 0.0241 Inf   2.138  0.1800
# pronoun - gap long        KLE    -0.0178 0.0265 Inf  -0.674  0.9844
# pronoun - gap island      KLE     0.0269 0.0228 Inf   1.184  0.8017
# pronoun - gap short       MLE    -0.0275 0.0230 Inf  -1.194  0.7956
# pronoun - gap long        MLE    -0.0265 0.0252 Inf  -1.050  0.8755
# pronoun - gap island      MLE    -0.0417 0.0217 Inf  -1.919  0.2876

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()
  
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment) %>%
  mutate(across(4:6, round, 2), across(7, round, 3))

dt %>% write_rds('tables/orc_spr_logrt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_spr_logrt_proficiency_posthoc.rds')

# plot ...

emmip(model3, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4 <- lmer(log(rt) ~ dependency * environment * group * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 + dependency + environment + group | item),
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model4)
toc()
beep()

# ... 75.33 sec elapsed

performance::check_convergence(model4) # convergence: okay
check_singularity(model4, tolerance = 1e-9) # singularity: warning
check_collinearity(model4) # independence: okay
check_heteroscedasticity(model4) # homogeneity of variance: okay
check_outliers(model4) # outliers: okay
check_autocorrelation(model4) # autocorrelation: warning
MuMIn::r.squaredGLMM(model4) # R2m = 0.07532366, R2c = 0.3664508
check_model(model4) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15)

# save ...

model4 %>% write_rds('models/src_spr_logrt_proficiency_model4.rds')
model4 <- read_rds('models/src_spr_logrt_proficiency_model4.rds')

# summary ...

summary(model4)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                  7.272e+00  2.197e-02  1.312e+02 330.974  < 2e-16 ***
# dependency2                                 -7.609e-02  1.398e-02  6.725e+01  -5.445 7.89e-07 ***
# environment2                                 1.874e-02  1.730e-02  3.696e+01   1.083   0.2857    
# environment3                                -9.478e-02  1.750e-02  5.763e+01  -5.414 1.25e-06 ***
# group2                                       2.005e-01  4.330e-02  1.292e+02   4.632 8.73e-06 ***
# proficiency                                  1.745e-03  2.184e-02  1.292e+02   0.080   0.9364    
# dependency2:environment2                    -1.323e-01  2.904e-02  3.481e+03  -4.556 5.38e-06 ***
# dependency2:environment3                    -1.237e-01  2.894e-02  3.485e+03  -4.275 1.96e-05 ***
# dependency2:group2                          -3.948e-02  2.673e-02  1.451e+02  -1.477   0.1418    
# environment2:group2                          5.013e-02  3.192e-02  1.722e+02   1.570   0.1182    
# environment3:group2                          2.239e-02  3.417e-02  1.438e+02   0.655   0.5134    
# dependency2:proficiency                     -1.822e-02  1.354e-02  1.472e+02  -1.345   0.1806    
# environment2:proficiency                     3.756e-03  1.621e-02  1.758e+02   0.232   0.8170    
# environment3:proficiency                     1.111e-02  1.730e-02  1.455e+02   0.642   0.5219    
# group2:proficiency                           2.743e-02  4.369e-02  1.293e+02   0.628   0.5311    
# dependency2:environment2:group2             -6.286e-02  5.800e-02  3.507e+03  -1.084   0.2786    
# dependency2:environment3:group2             -4.290e-02  5.788e-02  3.504e+03  -0.741   0.4586    
# dependency2:environment2:proficiency        -6.728e-02  2.947e-02  3.498e+03  -2.283   0.0225 *  
# dependency2:environment3:proficiency        -1.653e-02  2.931e-02  3.497e+03  -0.564   0.5729    
# dependency2:group2:proficiency               2.669e-03  2.710e-02  1.473e+02   0.098   0.9217    
# environment2:group2:proficiency             -3.375e-02  3.249e-02  1.769e+02  -1.039   0.3003    
# environment3:group2:proficiency             -1.591e-02  3.466e-02  1.463e+02  -0.459   0.6468    
# dependency2:environment2:group2:proficiency  3.826e-02  5.899e-02  3.526e+03   0.649   0.5166    
# dependency2:environment3:group2:proficiency -1.753e-02  5.868e-02  3.522e+03  -0.299   0.7652  

# table ...

dt <- summary(model4)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3')) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_proficiency_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                            | Estimate|   SE|      t|     p|
# |:----------------------------------------------------|--------:|----:|------:|-----:|
# |(Intercept)                                          |     7.27| 0.02| 330.97| 0.000| ***
# |Dependency[RP]                                       |    -0.08| 0.01|  -5.44| 0.000| ***
# |Environment[Long]                                    |     0.02| 0.02|   1.08| 0.286|
# |Environment[Island]                                  |    -0.09| 0.02|  -5.41| 0.000| ***
# |Group[KLE]                                           |     0.20| 0.04|   4.63| 0.000| ***
# |Proficiency                                          |     0.00| 0.02|   0.08| 0.936|
# |Dependency[RP]:Environment[Long]                     |    -0.13| 0.03|  -4.56| 0.000| ***
# |Dependency[RP]:Environment[Island]                   |    -0.12| 0.03|  -4.28| 0.000| ***
# |Dependency[RP]:Group[KLE]                            |    -0.04| 0.03|  -1.48| 0.142|
# |Environment[Long]:Group[KLE]                         |     0.05| 0.03|   1.57| 0.118|
# |Environment[Island]:Group[KLE]                       |     0.02| 0.03|   0.66| 0.513|
# |Dependency[RP]:Proficiency                           |    -0.02| 0.01|  -1.35| 0.181|
# |Environment[Long]:Proficiency                        |     0.00| 0.02|   0.23| 0.817|
# |Environment[Island]:Proficiency                      |     0.01| 0.02|   0.64| 0.522|
# |Group[KLE]:Proficiency                               |     0.03| 0.04|   0.63| 0.531|
# |Dependency[RP]:Envir.[Long]:Group[KLE]               |    -0.06| 0.06|  -1.08| 0.279|
# |Dependency[RP]:Envir.[Island]:Group[KLE]             |    -0.04| 0.06|  -0.74| 0.459|
# |Dependency[RP]:Envir.[Long]:Proficiency              |    -0.07| 0.03|  -2.28| 0.022| *
# |Dependency[RP]:Envir.[Island]:Proficiency            |    -0.02| 0.03|  -0.56| 0.573|
# |Dependency[RP]:Group[KLE]:Proficiency                |     0.00| 0.03|   0.10| 0.922|
# |Environment[Long]:Group[KLE]:Proficiency             |    -0.03| 0.03|  -1.04| 0.300|
# |Environment[Island]:Group[KLE]:Proficiency           |    -0.02| 0.03|  -0.46| 0.647|
# |Dependency[RP]:Envir.[Long]:Group[KLE]:Proficiency   |     0.04| 0.06|   0.65| 0.517|
# |Dependency[RP]:Envir.[Island]:Group[KLE]:Proficiency |    -0.02| 0.06|  -0.30| 0.765|


# joint tests ...

joint_tests(model4)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  29.716  <.0001 ***
# environment                                2 Inf  20.738  <.0001 ***
# group                                      1 Inf  21.484  <.0001 ***
# proficiency                                1 Inf   0.006  0.9363
# dependency:environment                     2 Inf  13.076  <.0001 ***
# dependency:group                           1 Inf   2.181  0.1397
# dependency:proficiency                     1 Inf   1.810  0.1785
# environment:group                          2 Inf   1.226  0.2933
# environment:proficiency                    2 Inf   0.206  0.8135
# group:proficiency                          1 Inf   0.394  0.5300
# dependency:environment:group               2 Inf   0.611  0.5427
# dependency:environment:proficiency         2 Inf   2.832  0.0589 .
# dependency:group:proficiency               1 Inf   0.010  0.9215
# environment:group:proficiency              2 Inf   0.541  0.5824
# dependency:environment:group:proficiency   2 Inf   0.471  0.6244

# post-hoc tests ...

pairwise <- model4 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       KLE    0.01184 0.0321 Inf   0.368  0.9994
# pronoun - gap long        KLE   -0.07458 0.0323 Inf  -2.306  0.1200
# pronoun - gap island      KLE    0.00408 0.0321 Inf   0.127  1.0000
# pronoun - gap short       MLE    0.00759 0.0293 Inf   0.259  0.9999
# pronoun - gap long        MLE   -0.04056 0.0292 Inf  -1.387  0.6606
# pronoun - gap island      MLE   -0.01770 0.0288 Inf  -0.615  0.9902

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment) %>%
  mutate(across(4:6, round, 2), across(7, round, 3))

table %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|     z|     p|
# |:-----|:-----------|:--------|--------:|----:|-----:|-----:|
# |KLE   |Short       |RP - Gap |     0.01| 0.03|  0.37| 0.999|
# |KLE   |Long        |RP - Gap |    -0.07| 0.03| -2.31| 0.120|
# |KLE   |Island      |RP - Gap |     0.00| 0.03|  0.13| 1.000|
# |MLE   |Short       |RP - Gap |     0.01| 0.03|  0.26| 1.000|
# |MLE   |Long        |RP - Gap |    -0.04| 0.03| -1.39| 0.661|
# |MLE   |Island      |RP - Gap |    -0.02| 0.03| -0.61| 0.990|

table %>% write_rds('tables/src_spr_logrt_proficiency_posthoc.rds')
csv <- table %>% csv_table() %>% write_csv('tables/src_spr_logrt_proficiency_posthoc.csv')
table <- read_rds('tables/src_spr_logrt_proficiency_posthoc.rds')

table <- table %>% my_table()

# plot ...

emmip(model4, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: raw RTs (combined critical region) ####
# ---------------------------------------------------------------------------- #

# calculate rts for combined critical region ...

spr_combo <- spr_rt %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  group_by(study, group, participant, item, dependency, environment, accuracy, proficiency) %>%
  summarise(rt = sum(rt)) %>%
  ungroup()

# replacements ...

spr_combo <- spr_combo %>%
  group_by(study, group, dependency, environment) %>%
  mutate(sd2 = mean(rt, na.rm = T) + (2 * (sd(rt, na.rm = T)))) %>%
  ungroup() %>%
  mutate(rt = case_when(rt > sd2 ~ sd2, TRUE ~ as.numeric(rt))) %>%
  select(-sd2)
  
# summarise ...

plot <- spr_combo %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  mutate(environment = as.factor(environment)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean_rt = mean(rt, na.rm=T),
            lcl = my.lcl(rt),
            ucl = my.ucl(rt)) %>%
  ungroup() %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean_rt, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl, col = dependency), 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1,
             position = position_dodge(width = .4),
             shape = 21,
             fill = 'white',
             show.legend = FALSE),
  theme_classic(),
  scale_y_continuous(name='RT (ms)'),
  scale_x_discrete(name='Environment'),
  scale_colour_manual(name='Dependency:', 
                      values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        legend.position = 'bottom',
        plot.margin = margin(2, 2, 0, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_wrap2(vars(group_label), axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_spr_rt_combo.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_spr_rt_combo.png', width=6.5, height=3, dpi=600)

# table ...

table <- plot %>%
  select(study, Group = group, Region = region2, Environment = environment, 
         Dependency = dependency, Mean = mean_rt, LCL = lcl, UCL = ucl) %>%
  arrange(study, Group, Region, Environment, Dependency) %>%
  filter(Region %in% c(1, 2, 3)) %>%
  mutate(Dependency = as.character(Dependency)) %>%
  mutate_if(is.character, str_replace_all, 
            c('gap' = 'Gap', 'pronoun' = 'RP'))

table %>% write_rds('tables/spr_rt_combo.rds')

# ---------------------------------------------------------------------------- #
# plot: RRTs
# ---------------------------------------------------------------------------- #

# summarise ...

plot <- spr_rrt %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  mutate(environment = as.factor(environment)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, region2, dependency, environment, condition) %>%
  summarise(mean_rrt = mean(rrt, na.rm=T),
            lcl = my.lcl(rrt),
            ucl = my.ucl(rrt),
            sd = sd(rrt, na.rm=T)) %>%
  ungroup() %>%
  filter(region2 %in% c(-2, -1, 0, 1, 2, 3, 4)) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island')) %>%
  group_by(region2) %>%
  mutate(width = 0.01 * n()) %>%
  ungroup()

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=region2, y=mean_rrt, group=dependency, col=dependency, shape=dependency, width=width*1.5),
  geom_hline(yintercept = 0),
  geom_vline(xintercept = 0),
  annotate('rect', xmin = 0.5, xmax = 3.5, ymin = min(plot$lcl) - .01, ymax = max(plot$ucl) + .01, alpha = 0, col = 'grey75', linetype = 2),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl, col = dependency), 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1,
             position = position_dodge(width = .4),
             shape = 21,
             fill = 'white',
             show.legend = FALSE),
  theme_classic(),
  ylab('RRT (ms)'),
  # scale_y_continuous(name="RRT (ms)", breaks = c(-100, 0, 100), limits=c(-164, 164)),
  scale_x_continuous(name='Region', limits=c(-2.25, 4.25), breaks = c(-2, -1, 0, 1, 2, 3, 4)),
  scale_colour_manual(name='Dependency:', values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        legend.position = 'bottom',
        plot.margin = margin(2, 2, 0, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(vars(group_label), vars(environment), axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_spr_rrt.png', width=6.5, height=4, dpi=600)

p2 + s
ggsave('plots/src_spr_rrt.png', width=6.5, height=4, dpi=600)

# table ...

table <- plot %>%
  select(study, Group = group, Region = region2, Environment = environment, 
         Dependency = dependency, Mean = mean_rrt, LCL = lcl, UCL = ucl) %>%
  arrange(study, Group, Region, Environment, Dependency) %>%
  filter(Region %in% c(1, 2, 3)) %>%
  mutate(Dependency = as.character(Dependency)) %>%
  mutate_if(is.character, str_replace_all, 
            c('gap' = 'Gap', 'pronoun' = 'RP'))

table %>% write_rds('tables/spr_rrt.rds')

# ---------------------------------------------------------------------------- #
# model: RRT (old) 
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_rrt %>%
  select(study, group, participant, item, region2, dependency, environment, rrt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'region2', 'dependency', 'environment'), factor) %>%
  mutate(rrt = as.numeric(rrt)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc region 1 model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(rrt ~ dependency * environment * group + 
                 (1 + dependency | participant) + 
                 (1 + dependency | item), 
               data = filter(md, study == 'ORC', region2 == '1'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

performance::check_convergence(model3) # TRUE
performance::check_singularity(model3) # FALSE

# no warnings

model3 %>% write_rds('models/orc_spr_rrt_region1.rds')
model3 <- read_rds('models/orc_spr_rrt_region1.rds')

# model summary ...

summary(model3)

# Parameter                        Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                      -23.6140     5.9071  122.3432  -3.998 0.000110 ***
# dependency2                      -19.6894     6.6481   43.6300  -2.962 0.004936 ** 
# environment2                       2.2730     6.2705 6045.6242   0.362 0.717002    
# environment3                       7.0803     6.2546 6078.0329   1.132 0.257675    
# group2                           -32.7735    12.6748  222.2722  -2.586 0.010356 *  
# group3                           -23.5970    12.5022  221.2128  -1.887 0.060413 .  
# dependency2:environment2          -0.5114    12.5624 6094.7375  -0.041 0.967530    
# dependency2:environment3           9.4557    12.5093 6077.0363   0.756 0.449741    
# dependency2:group2                 0.5883    14.6778  220.7517   0.040 0.968066    
# dependency2:group3                 1.2924    14.4482  218.2985   0.089 0.928803    
# environment2:group2                6.6305    15.1539 6066.8044   0.438 0.661732    
# environment3:group2                9.1589    15.1314 6068.7811   0.605 0.545009    
# environment2:group3                6.5280    14.9063 6086.2952   0.438 0.661449    
# environment3:group3                4.8520    14.8949 6076.8528   0.326 0.744627    
# dependency2:environment2:group2   29.2793    30.3019 6069.5502   0.966 0.333956    
# dependency2:environment3:group2   80.9886    30.2596 6069.3640   2.676 0.007460 ** 
# dependency2:environment2:group3   19.0361    29.7999 6087.8756   0.639 0.522978    
# dependency2:environment3:group3  109.8254    29.7740 6075.8345   3.689 0.000227 ***

# table ...

table <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/orc_spr_rrt_region1_model.rds')

# joint tests ...

joint_tests(model3)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf   8.771  0.0031 **
# environment                    2 Inf   0.668  0.5125
# group                          2 Inf   3.705  0.0246 *
# dependency:environment         2 Inf   0.402  0.6689
# dependency:group               2 Inf   0.004  0.9960
# environment:group              4 Inf   0.118  0.9761
# dependency:environment:group   4 Inf   4.248  0.0019 **

# post-hoc tests ...

pairwise <- model3 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS       3.27 15.2 Inf   0.215  1.0000
# pronoun - gap long        ENS     -13.34 15.2 Inf  -0.879  0.9856
# pronoun - gap island      ENS     -50.88 15.2 Inf  -3.346  0.0073 **
# pronoun - gap short       KLE     -32.90 17.5 Inf  -1.882  0.4218
# pronoun - gap long        KLE     -20.23 17.5 Inf  -1.155  0.9202
# pronoun - gap island      KLE      -6.06 17.4 Inf  -0.348  1.0000
# pronoun - gap short       MLE     -38.39 17.0 Inf  -2.259  0.1941
# pronoun - gap long        MLE     -35.97 17.0 Inf  -2.119  0.2658
# pronoun - gap island      MLE      17.29 16.9 Inf   1.021  0.9616

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/orc_spr_rrt_region1_posthoc.rds')

# plot ...

plot(emmeans(model3, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... orc region 2 model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model2 <- lmer(rrt ~ dependency * environment * group + 
                 (1 + dependency | participant) + 
                 (1 | item), 
               data = filter(md, study == 'ORC', region2 == '2'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

performance::check_convergence(model2) # TRUE
performance::check_singularity(model2) # FALSE

# no warnings

model2 %>% write_rds('models/orc_spr_rrt_region2.rds')
model2 <- read_rds('models/orc_spr_rrt_region2.rds')

# model summary ...

summary(model2)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -20.430      6.699  214.601  -3.050 0.002580 ** 
# dependency2                      -30.833      4.676  223.618  -6.594 3.05e-10 ***
# environment2                       5.475      4.766 6001.170   1.149 0.250707    
# environment3                      13.068      4.770 5999.094   2.740 0.006168 ** 
# group2                           -70.784     15.930  224.168  -4.443 1.39e-05 ***
# group3                           -74.942     15.719  223.356  -4.768 3.36e-06 ***
# dependency2:environment2         -29.484      9.545 5972.971  -3.089 0.002018 ** 
# dependency2:environment3         -30.517      9.540 5989.062  -3.199 0.001387 ** 
# dependency2:group2               -43.358     11.328  227.413  -3.827 0.000167 ***
# dependency2:group3               -55.229     11.113  221.609  -4.970 1.34e-06 ***
# environment2:group2               20.979     11.568 6006.191   1.814 0.069801 .  
# environment3:group2               14.725     11.581 6002.208   1.271 0.203612    
# environment2:group3               26.716     11.321 5999.161   2.360 0.018311 *  
# environment3:group3               19.901     11.332 5997.468   1.756 0.079098 .  
# dependency2:environment2:group2  -25.588     23.132 5984.735  -1.106 0.268691    
# dependency2:environment3:group2   30.897     23.164 5980.922   1.334 0.182312    
# dependency2:environment2:group3  -11.987     22.635 5991.119  -0.530 0.596425    
# dependency2:environment3:group3   28.280     22.652 5979.303   1.248 0.211921  

# table ...

table <- summary(model2)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/orc_spr_rrt_region2_model.rds')

# joint tests ...

joint_tests(model2)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  43.476  <.0001 ***
# environment                    2 Inf   3.784  0.0227 *
# group                          2 Inf  14.847  <.0001 ***
# dependency:environment         2 Inf   6.605  0.0014 **
# dependency:group               2 Inf  14.036  <.0001 ***
# environment:group              4 Inf   1.696  0.1478
# dependency:environment:group   4 Inf   1.683  0.1507

# post-hoc tests ...

pairwise <- model2 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS      24.43 11.5 Inf   2.124  0.2632
# pronoun - gap long        ENS       7.47 11.5 Inf   0.652  0.9984
# pronoun - gap island      ENS     -25.81 11.5 Inf  -2.244  0.2015
# pronoun - gap short       KLE     -20.70 13.2 Inf  -1.571  0.6673
# pronoun - gap long        KLE     -63.25 13.3 Inf  -4.767  <.0001 ***
# pronoun - gap island      KLE     -40.04 13.3 Inf  -3.013  0.0230 *
# pronoun - gap short       MLE     -36.23 12.7 Inf  -2.846  0.0391 *
# pronoun - gap long        MLE     -65.18 12.8 Inf  -5.109  <.0001 ***
# pronoun - gap island      MLE     -58.19 12.8 Inf  -4.562  <.0001 ***

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/orc_spr_rrt_region2_posthoc.rds')

# plot ...

plot(emmeans(model2, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... orc region 3 model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model3 <- lmer(rrt ~ dependency * environment * group + 
                 (1 | participant) + 
                 (1 + dependency | item), 
               data = filter(md, study == 'ORC', region2 == '3'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

performance::check_convergence(model3) # TRUE
performance::check_singularity(model3) # FALSE

# no warnings

model3 %>% write_rds('models/orc_spr_rrt_region3.rds')
model3 <- read_rds('models/orc_spr_rrt_region3.rds')

# model summary ...

summary(model3)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -18.866      6.927  168.264  -2.724 0.007141 ** 
# dependency2                      -12.139      4.267   28.319  -2.845 0.008160 ** 
# environment2                     -10.461      5.100 5867.156  -2.051 0.040308 *  
# environment3                       4.346      5.095 5920.657   0.853 0.393714    
# group2                           -70.629     14.828  225.036  -4.763 3.41e-06 ***
# group3                           -54.237     14.579  221.321  -3.720 0.000252 ***
# dependency2:environment2           5.805     10.221 5946.047   0.568 0.570046    
# dependency2:environment3          -1.259     10.184 5927.351  -0.124 0.901619    
# dependency2:group2               -33.245     10.153 5904.164  -3.274 0.001065 ** 
# dependency2:group3               -23.602      9.752 5900.099  -2.420 0.015545 *  
# environment2:group2              -10.672     12.426 5910.759  -0.859 0.390481    
# environment3:group2                5.490     12.438 5911.498   0.441 0.658934    
# environment2:group3               -6.970     11.979 5935.838  -0.582 0.560689    
# environment3:group3               -4.024     11.961 5924.436  -0.336 0.736536    
# dependency2:environment2:group2    2.861     24.832 5918.467   0.115 0.908293    
# dependency2:environment3:group2   24.469     24.873 5924.455   0.984 0.325267    
# dependency2:environment2:group3    9.628     23.933 5930.112   0.402 0.687492    
# dependency2:environment3:group3   22.402     23.887 5922.993   0.938 0.348368  

# table ...

table <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/orc_spr_rrt_region3_model.rds')

# joint tests ...

joint_tests(model3)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf   8.095  0.0044 **
# environment                    2 Inf   4.472  0.0114 *
# group                          2 Inf  13.040  <.0001 ***
# dependency:environment         2 Inf   0.272  0.7619
# dependency:group               2 Inf   5.984  0.0025 **
# environment:group              4 Inf   0.522  0.7195
# dependency:environment:group   4 Inf   0.374  0.8276

# post-hoc tests ...

pairwise <- model3 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS      11.89 11.3 Inf   1.052  0.9556
# pronoun - gap long        ENS      13.53 11.4 Inf   1.186  0.9109
# pronoun - gap island      ENS      -4.99 11.4 Inf  -0.436  0.9999
# pronoun - gap short       KLE     -30.47 13.5 Inf  -2.250  0.1999
# pronoun - gap long        KLE     -25.96 13.5 Inf  -1.928  0.3925
# pronoun - gap island      KLE     -22.88 13.5 Inf  -1.696  0.5713
# pronoun - gap short       MLE     -22.39 12.7 Inf  -1.764  0.5174
# pronoun - gap long        MLE     -11.12 12.5 Inf  -0.887  0.9855
# pronoun - gap island      MLE     -16.87 12.5 Inf  -1.354  0.8242

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/orc_spr_rrt_region3_posthoc.rds')

# plot ...

plot(emmeans(model3, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src region 1 model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(rrt ~ dependency * environment * group + 
                 (1 | participant) + 
                 (1 | item), 
               data = filter(md, study == 'SRC', region2 == '1'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

performance::check_convergence(model1) # TRUE
performance::check_singularity(model1) # FALSE

# no warnings

model1 %>% write_rds('models/src_spr_rrt_region1.rds')
model1 <- read_rds('models/src_spr_rrt_region1.rds')

# model summary ...

summary(model1)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -19.036      6.415   72.693  -2.967 0.004064 ** 
# dependency2                      -31.489      8.329 5238.298  -3.781 0.000158 ***
# environment2                      -7.970     10.233 5221.511  -0.779 0.436073    
# environment3                      16.310     10.194 5226.702   1.600 0.109684    
# group2                            66.614     15.370  170.905   4.334 2.49e-05 ***
# group3                            95.169     15.101  168.607   6.302 2.47e-09 ***
# dependency2:environment2         -43.947     20.479 5248.455  -2.146 0.031921 *  
# dependency2:environment3         -69.225     20.387 5233.107  -3.395 0.000690 ***
# dependency2:group2               -39.918     20.797 5250.079  -1.919 0.054983 .  
# dependency2:group3              -101.003     20.348 5225.232  -4.964 7.13e-07 ***
# environment2:group2              -14.900     25.593 5223.271  -0.582 0.560460    
# environment3:group2              -12.921     25.450 5229.237  -0.508 0.611691    
# environment2:group3               -1.184     24.967 5234.787  -0.047 0.962174    
# environment3:group3               39.976     24.884 5234.832   1.607 0.108219    
# dependency2:environment2:group2  -11.826     51.208 5249.740  -0.231 0.817376    
# dependency2:environment3:group2   36.681     50.908 5251.363   0.721 0.471226    
# dependency2:environment2:group3 -104.683     49.925 5220.800  -2.097 0.036059 *  
# dependency2:environment3:group3  -54.956     49.759 5223.470  -1.104 0.269452  

# table ...

table <- summary(model1)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/src_spr_rrt_region1_model.rds')

# joint tests ...

joint_tests(model1)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  14.294  0.0002 ***
# environment                    2 Inf   2.966  0.0515 *
# group                          2 Inf  20.585  <.0001 ***
# dependency:environment         2 Inf   5.891  0.0028 **
# dependency:group               2 Inf  12.644  <.0001 ***
# environment:group              4 Inf   1.518  0.1938
# dependency:environment:group   4 Inf   1.704  0.1459

# post-hoc tests ...

pairwise <- model1 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS      38.23 25.9 Inf   1.477  0.7419
# pronoun - gap long        ENS      33.12 25.9 Inf   1.277  0.8683
# pronoun - gap island      ENS     -24.90 25.6 Inf  -0.973  0.9731
# pronoun - gap short       KLE      -9.97 25.5 Inf  -0.391  1.0000
# pronoun - gap long        KLE     -26.91 25.1 Inf  -1.074  0.9498
# pronoun - gap island      KLE     -36.42 24.8 Inf  -1.470  0.7472
# pronoun - gap short       MLE      -9.56 23.9 Inf  -0.399  1.0000
# pronoun - gap long        MLE    -119.35 24.0 Inf  -4.970  <.0001 ***
# pronoun - gap island      MLE    -127.65 24.0 Inf  -5.312  <.0001 ***

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/src_spr_rrt_region1_posthoc.rds')

# plot ...

plot(emmeans(model1, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src region 2 model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(rrt ~ dependency * environment * group + 
                 (1 | participant) + 
                 (1 | item), 
               data = filter(md, study == 'SRC', region2 == '2'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

performance::check_convergence(model1) # TRUE
performance::check_singularity(model1) # FALSE

# no warnings

model1 %>% write_rds('models/src_spr_rrt_region2.rds')
model1 <- read_rds('models/src_spr_rrt_region2.rds')

# model summary ...

summary(model1)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -22.678      5.568   58.790  -4.073 0.000141 ***
# dependency2                      -47.524      7.044 5368.984  -6.747 1.67e-11 ***
# environment2                      33.092      8.627 5360.583   3.836 0.000127 ***
# environment3                      19.886      8.632 5363.161   2.304 0.021276 *  
# group2                            15.928     12.120  195.212   1.314 0.190325    
# group3                             1.923     11.889  192.194   0.162 0.871689    
# dependency2:environment2        -103.085     17.279 5385.132  -5.966 2.59e-09 ***
# dependency2:environment3         -49.438     17.265 5364.470  -2.864 0.004206 ** 
# dependency2:group2               -30.034     17.621 5382.696  -1.704 0.088361 .  
# dependency2:group3               -59.273     17.216 5358.836  -3.443 0.000580 ***
# environment2:group2               30.443     21.552 5363.114   1.413 0.157852    
# environment3:group2              -44.678     21.591 5366.241  -2.069 0.038568 *  
# environment2:group3               90.176     21.100 5370.453   4.274 1.96e-05 ***
# environment3:group3              -10.385     21.097 5373.840  -0.492 0.622565    
# dependency2:environment2:group2  -63.918     43.130 5379.145  -1.482 0.138400    
# dependency2:environment3:group2  108.983     43.212 5381.689   2.522 0.011696 *  
# dependency2:environment2:group3 -145.167     42.185 5358.843  -3.441 0.000583 ***
# dependency2:environment3:group3   54.560     42.169 5356.452   1.294 0.195779 

# joint tests ...

joint_tests(model1)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  45.516  <.0001 ***
# environment                    2 Inf   7.453  0.0006 ***
# group                          2 Inf   1.057  0.3474
# dependency:environment         2 Inf  17.806  <.0001 ***
# dependency:group               2 Inf   5.935  0.0026 **
# environment:group              4 Inf   8.002  <.0001 ***
# dependency:environment:group   4 Inf   7.321  <.0001 ***

# table ...

table <- summary(model1)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/src_spr_rrt_region2_model.rds')

# post-hoc tests ...

pairwise <- model1 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS      28.03 21.9 Inf   1.281  0.8662
# pronoun - gap long        ENS      -5.36 21.9 Inf  -0.245  1.0000
# pronoun - gap island      ENS     -75.93 22.0 Inf  -3.458  0.0049 **
# pronoun - gap short       KLE     -17.03 21.3 Inf  -0.798  0.9931
# pronoun - gap long        KLE    -114.34 21.2 Inf  -5.404  <.0001 ***
# pronoun - gap island      KLE     -12.00 21.2 Inf  -0.565  0.9995
# pronoun - gap short       MLE      -1.05 20.3 Inf  -0.051  1.0000
# pronoun - gap long        MLE    -179.60 20.2 Inf  -8.876  <.0001 ***
# pronoun - gap island      MLE     -50.44 20.1 Inf  -2.504  0.1052

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/src_spr_rrt_region2_posthoc.rds')

# plot ...

plot(emmeans(model1, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src region 3 model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- lmer(rrt ~ dependency * environment * group + 
                 (1 + environment | participant) + 
                 (1 + dependency | item), 
               data = filter(md, study == 'SRC', region2 == '3'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model7)
toc()
beep()

performance::check_convergence(model7) # TRUE
performance::check_singularity(model7) # FALSE

# no warnings

model7 %>% write_rds('models/src_spr_rrt_region3.rds')
model7 <- read_rds('models/src_spr_rrt_region3.rds')

# model summary ...

summary(model7)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -19.329      6.167   63.495  -3.134  0.00261 ** 
# dependency2                      -16.330      8.587   26.696  -1.902  0.06805 .  
# environment2                      17.121      9.452  160.508   1.811  0.07196 .  
# environment3                      58.738     13.554  175.769   4.334 2.46e-05 ***
# group2                            17.491     13.148  178.574   1.330  0.18513    
# group3                            40.031     12.935  176.930   3.095  0.00229 ** 
# dependency2:environment2         -37.049     17.186 4906.509  -2.156  0.03115 *  
# dependency2:environment3         -74.494     17.149 4887.720  -4.344 1.43e-05 ***
# dependency2:group2                 1.175     17.470 4897.623   0.067  0.94639    
# dependency2:group3               -14.493     17.130 4875.818  -0.846  0.39756    
# environment2:group2               -2.590     23.506  160.374  -0.110  0.91239    
# environment3:group2              -89.464     33.766  176.991  -2.650  0.00879 ** 
# environment2:group3               19.798     23.168  159.572   0.855  0.39407    
# environment3:group3             -173.519     33.209  175.417  -5.225 4.90e-07 ***
# dependency2:environment2:group2  -95.823     42.731 4898.929  -2.242  0.02498 *  
# dependency2:environment3:group2  -73.551     42.831 4899.612  -1.717  0.08600 .  
# dependency2:environment2:group3    9.369     42.089 4900.246   0.223  0.82386    
# dependency2:environment3:group3   58.001     41.956 4891.833   1.382  0.16691 

# table ...

table <- summary(model7)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/src_spr_rrt_region3_model.rds')

# joint tests ...

joint_tests(model7)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf   3.617  0.0572 *
# environment                    2 Inf  10.264  <.0001 ***
# group                          2 Inf   4.850  0.0078 **
# dependency:environment         2 Inf   9.436  0.0001 ***
# dependency:group               2 Inf   0.538  0.5838
# environment:group              4 Inf   7.483  <.0001 ***
# dependency:environment:group   4 Inf   3.210  0.0121 *

# post-hoc tests ...

pairwise <- model7 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS      13.96 22.3 Inf   0.625  0.9989
# pronoun - gap long        ENS       5.73 22.2 Inf   0.258  1.0000
# pronoun - gap island      ENS     -55.35 22.3 Inf  -2.486  0.1101
# pronoun - gap short       KLE      71.59 21.6 Inf   3.316  0.0082 **
# pronoun - gap long        KLE     -32.46 21.6 Inf  -1.502  0.7212
# pronoun - gap island      KLE     -71.27 21.8 Inf  -3.274  0.0095 **
# pronoun - gap short       MLE     -22.99 21.0 Inf  -1.094  0.9427
# pronoun - gap long        MLE     -21.86 20.9 Inf  -1.045  0.9566
# pronoun - gap island      MLE     -34.30 20.6 Inf  -1.668  0.5912

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/src_spr_rrt_region3_posthoc.rds')

# plot ...

plot(emmeans(model7, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: RRT collapsed 
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_trim %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  group_by(study, group, participant, item, dependency, environment) %>%
  summarise(rrt = sum(rrt)) %>%
  ungroup() %>%
  select(study, group, participant, item, dependency, environment, rrt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(rrt = as.numeric(rrt)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean_rrt = mean(rrt, na.rm=T),
            lcl = my.lcl(rrt),
            ucl = my.ucl(rrt),
            sd = sd(rrt, na.rm=T)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

s <- list(
  aes(x=environment, y=mean_rrt, group=dependency, col=dependency, shape=dependency),
  #geom_hline(yintercept = 0),
  geom_line(position = position_dodge(width = .4)),
  geom_point(size = 3, alpha = .5, position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), linetype = 1, width = .5, position = position_dodge(width = .4)),
  theme_classic(),
  scale_y_continuous(name="RRT (ms)"),
  scale_x_discrete(name="Environment"),
  scale_colour_manual(name="Dependency", values=c('#648fff', '#ffb000'), labels=c('Gap', 'RP')),
  scale_shape_manual(name="Dependency", values=c(16, 15), labels=c('Gap', 'RP')),
  theme(
    text = element_text(size = 12),
    #legend.position = 'right',
    plot.margin = margin(2, 2, 2, 2),
    #legend.margin = margin(r = -.4, unit = 'cm')
  ),
  facet_manual(~group_label, design = 'ABC', axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_spr_rrt_collapsed.png', width=6.5, height=2, dpi=600)

p2 + s
ggsave('plots/src_spr_rrt_collapsed.png', width=6.5, height=2, dpi=600)

# table ...

table <- plot %>%
  select(study, Group = group, Environment = environment, 
         Dependency = dependency, Mean = mean_rrt, LCL = lcl, UCL = ucl) %>%
  arrange(study, Group, Environment, Dependency) %>%
  mutate(Dependency = as.character(Dependency)) %>%
  mutate_if(is.character, str_replace_all, 
            c('gap' = 'Gap', 'pronoun' = 'RP'))

table %>% write_rds('tables/spr_rrt_collapsed.rds')

# ---------------------------------------------------------------------------- #
# ... orc model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- lmer(rrt ~ dependency * environment * group + 
                 (1 + dependency | participant) + 
                 (1 + dependency + group | item), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model7)
toc()
beep()

performance::check_convergence(model7) # TRUE
performance::check_singularity(model7) # FALSE (lme4 says it's singular)

# singular fit?

summary(model7)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -56.799     15.814  215.027  -3.592 0.000407 ***
# dependency2                      -59.602     10.545  104.671  -5.652 1.38e-07 ***
# environment2                      -1.482     10.197 6212.214  -0.145 0.884439    
# environment3                      23.460     10.192 6221.272   2.302 0.021382 *  
# group2                          -160.586     36.089  215.570  -4.450 1.38e-05 ***
# group3                          -142.123     36.318  219.409  -3.913 0.000121 ***
# dependency2:environment2         -21.206     20.416 6217.918  -1.039 0.298984    
# dependency2:environment3         -21.857     20.330 6218.653  -1.075 0.282372    
# dependency2:group2               -66.360     24.436  225.218  -2.716 0.007128 ** 
# dependency2:group3               -76.125     24.103  224.141  -3.158 0.001805 ** 
# environment2:group2               19.429     24.569 6217.950   0.791 0.429100    
# environment3:group2               27.846     24.558 6220.171   1.134 0.256869    
# environment2:group3               25.660     24.295 6187.520   1.056 0.290934    
# environment3:group3               17.651     24.294 6186.484   0.727 0.467518    
# dependency2:environment2:group2   13.590     49.163 5816.501   0.276 0.782224    
# dependency2:environment3:group2  128.489     49.080 6209.794   2.618 0.008868 ** 
# dependency2:environment2:group3   13.038     48.395 6226.292   0.269 0.787630    
# dependency2:environment3:group3  150.884     48.382 6217.324   3.119 0.001825 ** 

joint_tests(model7)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  31.949  <.0001 ***
# environment                    2 Inf   3.767  0.0231 *
# group                          2 Inf  12.495  <.0001 ***
# dependency:environment         2 Inf   0.746  0.4740
# dependency:group               2 Inf   6.104  0.0022 **
# environment:group              4 Inf   0.502  0.7346
# dependency:environment:group   4 Inf   3.567  0.0065 **

pairwise <- model7 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS      36.25 24.6 Inf   1.471  0.7405
# pronoun - gap long        ENS       6.17 24.7 Inf   0.250  1.0000
# pronoun - gap island      ENS     -78.73 24.7 Inf  -3.190  0.0127 *
# pronoun - gap short       KLE     -77.47 28.6 Inf  -2.713  0.0581 *
# pronoun - gap long        KLE     -93.96 28.5 Inf  -3.296  0.0088 **
# pronoun - gap island      KLE     -63.96 28.5 Inf  -2.245  0.2004
# pronoun - gap short       MLE     -94.52 27.7 Inf  -3.407  0.0059 **
# pronoun - gap long        MLE    -111.56 27.8 Inf  -4.019  0.0005 ***
# pronoun - gap island      MLE     -58.62 27.8 Inf  -2.110  0.2707

# ---------------------------------------------------------------------------- #
# ... src model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- lmer(rrt ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 | item), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model7)
toc()
beep()

performance::check_convergence(model7) # TRUE
performance::check_singularity(model7) # FALSE

# no warnings

summary(model7)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                      -60.387     10.321   58.552  -5.851 2.33e-07 ***
# dependency2                      -93.059     14.345  192.980  -6.487 7.18e-10 ***
# environment2                      40.231     18.750  180.460   2.146 0.033234 *  
# environment3                      90.487     21.446  190.519   4.219 3.79e-05 ***
# group2                            91.813     24.718  193.699   3.714 0.000266 ***
# group3                           122.018     24.326  192.526   5.016 1.19e-06 ***
# dependency2:environment2        -175.519     33.876 4982.902  -5.181 2.29e-07 ***
# dependency2:environment3        -183.349     33.854 5010.411  -5.416 6.38e-08 ***
# dependency2:group2               -65.742     35.758  194.532  -1.839 0.067506 .  
# dependency2:group3              -172.192     35.134  192.487  -4.901 2.02e-06 ***
# environment2:group2                9.517     46.720  181.708   0.204 0.838813    
# environment3:group2             -135.299     53.410  191.585  -2.533 0.012104 *  
# environment2:group3               97.787     45.945  180.519   2.128 0.034664 *  
# environment3:group3             -131.158     52.539  190.397  -2.496 0.013395 *  
# dependency2:environment2:group2 -167.774     84.425 5014.290  -1.987 0.046948 *  
# dependency2:environment3:group2   70.777     84.408 5011.938   0.839 0.401781    
# dependency2:environment2:group3 -233.034     82.950 5005.825  -2.809 0.004984 ** 
# dependency2:environment3:group3   59.179     82.896 5009.409   0.714 0.475324  

joint_tests(model7)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  42.082  <.0001 ***
# environment                    2 Inf  10.156  <.0001 ***
# group                          2 Inf  13.386  <.0001 ***
# dependency:environment         2 Inf  18.731  <.0001 ***
# dependency:group               2 Inf  12.387  <.0001 ***
# environment:group              4 Inf   3.824  0.0041 **
# dependency:environment:group   4 Inf   3.805  0.0043 **

pairwise <- model7 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       ENS       75.8 43.5 Inf   1.741  0.5356
# pronoun - gap long        ENS       33.9 43.5 Inf   0.779  0.9942
# pronoun - gap island      ENS     -150.9 43.3 Inf  -3.481  0.0045 **
# pronoun - gap short       KLE       42.4 42.1 Inf   1.007  0.9664
# pronoun - gap long        KLE     -167.3 41.9 Inf  -3.990  0.0006 ***
# pronoun - gap island      KLE     -113.5 42.0 Inf  -2.703  0.0602 .
# pronoun - gap short       MLE      -38.5 40.5 Inf  -0.949  0.9770
# pronoun - gap long        MLE     -313.4 40.4 Inf  -7.748  <.0001 ***
# pronoun - gap island      MLE     -205.9 40.4 Inf  -5.091  <.0001 ***

# ---------------------------------------------------------------------------- #
# model: RRT proficiency 
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_trim %>%
  filter(group %in% c('KLE', 'MLE'), region2 %in% c(1, 2, 3)) %>%
  filter(is.na(proficiency) == FALSE) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  mutate_at(c('environment', 'dependency', 'group'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  group_by(study, group, group_label, participant, proficiency, item, dependency, environment, condition) %>%
  summarise(rrt = sum(rrt, na.rm=T)) %>%
  ungroup() %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant')) %>%
  filter(is.na(proficiency) == FALSE)

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group, group_label, participant, proficiency, dependency, environment) %>%
  summarise(rrt = mean(rrt, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = rrt, col = dependency, fill = dependency),
  geom_hline(yintercept = 0),
  geom_vline(xintercept = 0),
  geom_smooth(aes(linetype = dependency), method=lm, alpha = .75, linewidth = .5),
  #geom_point(shape = 1, alpha = .5),
  theme_classic(),
  scale_x_continuous(name='Proficiency (z-score)'),
  scale_y_continuous(name="RRT (ms)"),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                      values=c('#c1d2ff', '#ffdf99'), 
                      labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                    values=c('41', 'solid'), 
                    labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 12, hjust = .5), 
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(-750, 750))
ggsave('plots/orc_spr_rrt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(-500, 500))
ggsave('plots/src_spr_rrt_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- lmer(rrt ~ dependency * environment * group * proficiency + 
                 (1 + dependency | participant) + 
                 (1 + environment | item),
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model7)
toc()
beep()

performance::check_convergence(model7) # TRUE
performance::check_singularity(model7) # FALSE

# 4.25 sec elapsed
# no warnings

# model summary ...

summary(model7)

# Parameter                                   Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                 -108.414     22.351  137.068  -4.850 3.30e-06 ***
# dependency2                                  -82.700     15.040  133.775  -5.499 1.87e-07 ***
# environment2                                   3.801     16.228   25.156   0.234  0.81673    
# environment3                                  29.918     16.686   26.505   1.793  0.08439 .  
# group2                                        21.163     43.573  133.132   0.486  0.62798    
# proficiency                                   62.438     21.881  132.967   2.854  0.00502 ** 
# dependency2:environment2                     -15.666     29.302 3760.871  -0.535  0.59293    
# dependency2:environment3                      26.356     29.279 3743.120   0.900  0.36808    
# dependency2:group2                           -10.935     30.045  133.302  -0.364  0.71646    
# environment2:group2                            5.332     29.356 3732.769   0.182  0.85589    
# environment3:group2                          -11.954     29.346 3724.649  -0.407  0.68378    
# dependency2:proficiency                       -5.305     15.128  134.083  -0.351  0.72637    
# environment2:proficiency                       7.241     14.678 3758.913   0.493  0.62179    
# environment3:proficiency                      -2.359     14.685 3760.797  -0.161  0.87239    
# group2:proficiency                          -124.298     43.763  132.980  -2.840  0.00522 ** 
# dependency2:environment2:group2               -1.567     58.493 3751.321  -0.027  0.97862    
# dependency2:environment3:group2               24.302     58.464 3750.324   0.416  0.67767    
# dependency2:environment2:proficiency          -3.926     29.430 3712.911  -0.133  0.89388    
# dependency2:environment3:proficiency          -3.322     29.467 3738.688  -0.113  0.91024    
# dependency2:group2:proficiency               -39.748     30.201  133.405  -1.316  0.19039    
# environment2:group2:proficiency              -57.621     29.369 3762.761  -1.962  0.04984 *  
# environment3:group2:proficiency                3.943     29.384 3764.714   0.134  0.89325    
# dependency2:environment2:group2:proficiency   60.610     58.739 3768.682   1.032  0.30221    
# dependency2:environment3:group2:proficiency   11.399     58.795 3772.372   0.194  0.84628 

# table ...

table <- summary(model7)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'proficiency' = 'Proficiency',
              'group2' = 'Group[MLE]')) %>%
  mutate_if(is.character, str_replace_all,
            c(':Environment\\[Long\\]:' = ':Envir.\\[Long\\]:',
              ':Environment\\[Island\\]:' = ':Envir.\\[Island\\]:'))

table %>% write_rds('tables/orc_spr_rrt_proficiency_model.rds')

# joint tests ...

joint_tests(model7)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  30.228  <.0001 ***
# environment                                2 Inf   1.806  0.1643
# group                                      1 Inf   0.241  0.6237
# proficiency                                1 Inf   8.142  0.0043 **
# dependency:environment                     2 Inf   1.051  0.3497
# dependency:group                           1 Inf   0.131  0.7176
# dependency:proficiency                     1 Inf   0.123  0.7258
# environment:group                          2 Inf   0.184  0.8322
# environment:proficiency                    2 Inf   0.232  0.7926
# group:proficiency                          1 Inf   8.067  0.0045 **
# dependency:environment:group               2 Inf   0.123  0.8839
# dependency:environment:proficiency         2 Inf   0.010  0.9897
# dependency:group:proficiency               1 Inf   1.732  0.1881
# environment:group:proficiency              2 Inf   2.755  0.0636 .
# dependency:environment:group:proficiency   2 Inf   0.602  0.5480

# post-hoc tests ...

pairwise <- model7 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       KLE      28.99 32.9 Inf   0.880  0.9401
# pronoun - gap long        KLE      -5.24 32.9 Inf  -0.159  1.0000
# pronoun - gap island      KLE      19.96 33.0 Inf   0.606  0.9905
# pronoun - gap short       MLE     -34.76 31.3 Inf  -1.110  0.8398
# pronoun - gap long        MLE      -8.39 31.4 Inf  -0.267  0.9999
# pronoun - gap island      MLE     -32.39 31.4 Inf  -1.032  0.8802

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/orc_spr_rrt_proficiency_posthoc.rds')

# plot ...

emmip(model7, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) 
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model6 <- lmer(rrt ~ dependency * environment * group * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 | item),
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model6)
toc()
beep()

performance::check_convergence(model6) # TRUE
performance::check_singularity(model6) # FALSE

# no warnings

# model summary ...

summary(model6)

# Parameter                                   Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                  -21.456     14.229   53.635  -1.508  0.13744    
# dependency2                                 -131.789     20.250  130.533  -6.508 1.49e-09 ***
# environment2                                  52.004     25.956  127.524   2.003  0.04724 *  
# environment3                                  38.883     30.893  128.716   1.259  0.21045    
# group2                                        41.284     26.686  129.463   1.547  0.12430    
# proficiency                                   19.139     13.447  128.816   1.423  0.15707    
# dependency2:environment2                    -234.175     47.030 3410.215  -4.979 6.70e-07 ***
# dependency2:environment3                    -161.088     47.019 3426.714  -3.426  0.00062 ***
# dependency2:group2                          -101.473     40.502  130.549  -2.505  0.01346 *  
# environment2:group2                           69.017     51.918  127.620   1.329  0.18611    
# environment3:group2                            5.006     61.793  128.787   0.081  0.93555    
# dependency2:proficiency                        3.415     20.380  129.321   0.168  0.86719    
# environment2:proficiency                     -44.511     26.141  126.662  -1.703  0.09107 .  
# environment3:proficiency                      -4.935     31.123  127.952  -0.159  0.87427    
# group2:proficiency                            17.687     26.895  128.816   0.658  0.51194    
# dependency2:environment2:group2              -92.574     93.971 3429.396  -0.985  0.32463    
# dependency2:environment3:group2              -56.343     94.062 3426.972  -0.599  0.54921    
# dependency2:environment2:proficiency         -55.337     47.273 3415.755  -1.171  0.24184    
# dependency2:environment3:proficiency         -86.996     47.292 3416.293  -1.840  0.06592 .  
# dependency2:group2:proficiency               -10.401     40.769  129.439  -0.255  0.79904    
# environment2:group2:proficiency              -57.880     52.329  127.500  -1.106  0.27077    
# environment3:group2:proficiency              -78.208     62.289  128.473  -1.256  0.21155    
# dependency2:environment2:group2:proficiency   69.054     94.559 3419.710   0.730  0.46527    
# dependency2:environment3:group2:proficiency   10.491     94.613 3422.991   0.111  0.91172 

# table ...

table <- summary(model6)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'proficiency' = 'Proficiency',
              'group2' = 'Group[MLE]')) %>%
  mutate_if(is.character, str_replace_all,
            c(':Environment\\[Long\\]:' = ':Envir.\\[Long\\]:',
              ':Environment\\[Island\\]:' = ':Envir.\\[Island\\]:'))

table %>% write_rds('tables/src_spr_rrt_proficiency_model.rds')

# joint tests ...

joint_tests(model6)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  42.355  <.0001 ***
# environment                                2 Inf   2.567  0.0768 .
# group                                      1 Inf   2.396  0.1217
# proficiency                                1 Inf   2.026  0.1547
# dependency:environment                     2 Inf  12.985  <.0001 ***
# dependency:group                           1 Inf   6.279  0.0122 *
# dependency:proficiency                     1 Inf   0.028  0.8669
# environment:group                          2 Inf   0.884  0.4133
# environment:proficiency                    2 Inf   1.450  0.2346
# group:proficiency                          1 Inf   0.433  0.5108
# dependency:environment:group               2 Inf   0.492  0.6114
# dependency:environment:proficiency         2 Inf   1.734  0.1766
# dependency:group:proficiency               1 Inf   0.065  0.7986
# environment:group:proficiency              2 Inf   1.274  0.2796
# dependency:environment:group:proficiency   2 Inf   0.310  0.7335

# post-hoc tests ...

pairwise <- model6 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df z.ratio p.value
# pronoun - gap short       KLE       69.3 50.9 Inf   1.362  0.6804
# pronoun - gap long        KLE      -20.5 50.8 Inf  -0.405  0.9990
# pronoun - gap island      KLE      -22.9 50.9 Inf  -0.451  0.9982
# pronoun - gap short       MLE       32.4 45.4 Inf   0.714  0.9790
# pronoun - gap long        MLE       11.6 45.4 Inf   0.256  0.9999
# pronoun - gap island      MLE      -49.3 45.4 Inf  -1.088  0.8563

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/src_spr_rrt_proficiency_posthoc.rds')

# plot ...

emmip(model6, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: accuracy overall ####
# ---------------------------------------------------------------------------- #

plot <- spr_clean %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant) %>%
  summarise(acc_rate = mean(accuracy)) %>%
  ungroup()

p1 <- ggplot(data = filter(plot, study == 'ORC'))
p2 <- ggplot(data = filter(plot, study == 'SRC'))

s = list(
  aes(x = group_label, y = acc_rate * 100, fill = group, label = participant),
  geom_hline(yintercept = 50),
  geom_violin(),
  geom_boxplot(width = .1, fill='white'),
  theme_classic(),
  scale_x_discrete(name="Group"),
  scale_y_continuous(name="Accuracy (%)", 
                     limits=c(0, 100)),
  theme(text = element_text(size = 12), 
        plot.title = element_text(size = 12, hjust = .5), 
        legend.position = "hide")
)

p1 + s

p2 + s

# ---------------------------------------------------------------------------- #
# plot: accuracy critical ####
# ---------------------------------------------------------------------------- #

# summarize for plotting by group ...

plot <- spr_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(accuracy = as.logical(accuracy)) %>%
  group_by(study, group, group_label, participant, item, dependency, environment, accuracy) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean = mean(accuracy, na.rm=T) * 100,
            lcl = my.lcl(accuracy) * 100,
            ucl = my.ucl(accuracy) * 100) %>%
  ungroup()

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# create plot ...

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(50, 100)),
  scale_x_discrete(name="Environment", 
                   limits = c("short", "long", "island"), 
                   labels = c("Short", "Long", "Island")),
  scale_y_continuous(name="Accuracy (%)", 
                     breaks = c(0, 25, 50, 75, 100)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    strip.text = element_markdown(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin=margin(0, 0, 0, -5)),
  facet_manual(~group_label,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p1 + s
ggsave("plots/orc_spr_accuracy.png", width=6.5, height=2, dpi=600)

p2 + s
ggsave("plots/src_spr_accuracy.png", width=6.5, height=2, dpi=600)

# table ...

dt <- plot %>%
  select(-group_label) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  select(study, group, environment, 
         dependency, mean, lcl, ucl) %>%
  arrange(study, group, environment, dependency) %>%
  my_table()

dt %>% write_rds('tables/spr_accuracy.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_accuracy.csv')
dt <- read_rds('tables/spr_accuracy.rds')

# ---------------------------------------------------------------------------- #
# plot: accuracy iqr critical ####
# ---------------------------------------------------------------------------- #

# summarize for plotting by group ...

plot <- spr_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(accuracy = as.logical(accuracy)) %>%
  group_by(study, group, group_label, participant, dependency, environment) %>%
  summarise(accuracy = mean(accuracy)) %>%
  ungroup() %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(iqr = IQR(accuracy)) %>%
  ungroup()

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# create plot ...

s <- list(
  aes(x=environment, y=iqr, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  #coord_cartesian(ylim = c(50, 100)),
  scale_x_discrete(name="Environment", 
                   limits = c("short", "long", "island"), 
                   labels = c("Short", "Long", "Island")),
  scale_y_continuous(name="Accuracy (%)"),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    strip.text = element_markdown(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin=margin(0, 0, 0, -5)),
  facet_manual(~group_label,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p1 + s
ggsave("plots/orc_spr_accuracy_iqr.png", width=6.5, height=2, dpi=600)

p2 + s
ggsave("plots/src_spr_accuracy_iqr.png", width=6.5, height=2, dpi=600)

# table ...

dt <- plot %>%
  select(-group_label) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  select(study, group, environment, 
         dependency, mean, lcl, ucl) %>%
  arrange(study, group, environment, dependency) %>%
  my_table()

dt %>% write_rds('tables/spr_accuracy_iqr.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_accuracy_iqr.csv')
dt <- read_rds('tables/spr_accuracy_iqr.rds')

# ---------------------------------------------------------------------------- #
# plot: accuracy variance ####
# ---------------------------------------------------------------------------- #

# summarize for plotting by group ...

spr_var <- spr_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(accuracy = as.logical(accuracy)) %>%
  group_by(study, group, group_label, participant, dependency, environment, proficiency) %>%
  summarise(var = var(accuracy)) %>%
  ungroup()

plot <- spr_var %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean = mean(var),
            lcl = my.lcl(var),
            ucl = my.ucl(var)) %>%
  ungroup()

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# create plot ...

s <- list( 
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                width = 0.2, 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(0, .2)),
  scale_x_discrete(name="Environment", 
                   limits = c("short", "long", "island"), 
                   labels = c("Short", "Long", "Island")),
  scale_y_continuous(name="Variance", breaks = c(0, .1, .2)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    strip.text = element_markdown(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin=margin(0, 0, 0, -5)),
  facet_manual(~group_label,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p1 + s 
ggsave("plots/orc_spr_accuracy_var.png", width=6.5, height=2, dpi=600)

p2 + s 
ggsave("plots/src_spr_accuracy_var.png", width=6.5, height=2, dpi=600)

# table ...

dt <- plot %>%
  select(-group_label) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  select(study, group, environment, 
         dependency, mean, lcl, ucl) %>%
  arrange(study, group, environment, dependency) %>%
  my_table()

dt %>% write_rds('tables/spr_accuracy_var.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_accuracy_var.csv')
dt <- read_rds('tables/spr_accuracy_var.rds')

# ---------------------------------------------------------------------------- #
# model: accuracy variance ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_var %>%
  select(study, group, participant, dependency, environment, var, proficiency) %>%
  mutate_at(c('study', 'group', 'participant', 'dependency', 'environment'), factor) %>%
  mutate(var = as.numeric(var)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(var ~ dependency * environment * group + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 3 sec
# no warnings

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model

model <- model1

# save ...

model %>% write_rds('models/orc_spr_var_model.rds')
model <- read_rds('models/orc_spr_var_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_spr_var_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_var_model.csv')
dt <- read_rds('tables/orc_spr_var_model.rds')
csv <- read_csv('tables/orc_spr_var_model.csv')

dt %>% rename("Est." = 2) %>% kbl(format = 'markdown')

# |Parameter                                |Est.  |$SE$ |$t$   |$p$    |
# |:----------------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                              |0.06  |0.00 |17.22 |< .001 | ***
# |Dependency[RP]                           |−0.04 |0.00 |−8.82 |< .001 | ***
# |Environment[Long]                        |0.00  |0.00 |0.43  |0.667  |
# |Environment[Island]                      |0.00  |0.00 |0.48  |0.630  |
# |Group[KLE]                               |0.06  |0.01 |7.61  |< .001 | ***
# |Group[MLE]                               |0.02  |0.01 |2.96  |0.003  | **
# |Dependency[RP]:Environment[Long]         |0.00  |0.01 |−0.23 |0.822  |
# |Dependency[RP]:Environment[Island]       |−0.01 |0.01 |−0.86 |0.392  |
# |Dependency[RP]:Group[KLE]                |−0.05 |0.01 |−4.96 |< .001 | ***
# |Dependency[RP]:Group[MLE]                |−0.01 |0.01 |−0.78 |0.433  |
# |Environment[Long]:Group[KLE]             |0.00  |0.01 |0.00  |0.998  |
# |Environment[Island]:Group[KLE]           |0.00  |0.01 |−0.10 |0.919  |
# |Environment[Long]:Group[MLE]             |0.00  |0.01 |−0.42 |0.671  |
# |Environment[Island]:Group[MLE]           |0.00  |0.01 |0.24  |0.807  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |0.02  |0.02 |0.96  |0.338  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.01  |0.02 |0.47  |0.641  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |0.00  |0.02 |−0.01 |0.991  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |−0.03 |0.02 |−1.30 |0.193  |

# joint tests ...

joint_tests(model)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 224  77.707  <.0001 ***
# environment                    2 224   0.139  0.8704
# group                          2 224  29.024  <.0001 ***
# dependency:environment         2 448   0.394  0.6744
# dependency:group               2 224  13.424  <.0001 ***
# environment:group              4 224   0.167  0.9548
# dependency:environment:group   4 448   0.970  0.4235

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df t.ratio p.value
# pronoun - gap short       ENS    -0.0135 0.0105 667  -1.291  0.8598
# pronoun - gap long        ENS    -0.0225 0.0105 667  -2.150  0.2522
# pronoun - gap island      ENS    -0.0153 0.0105 667  -1.463  0.7511
# pronoun - gap short       KLE    -0.0726 0.0121 667  -5.980  <.0001 ***
# pronoun - gap long        KLE    -0.0605 0.0121 667  -4.984  <.0001 ***
# pronoun - gap island      KLE    -0.0641 0.0121 667  -5.283  <.0001 ***
# pronoun - gap short       MLE    -0.0116 0.0119 667  -0.977  0.9717
# pronoun - gap long        MLE    -0.0209 0.0119 667  -1.757  0.5231
# pronoun - gap island      MLE    -0.0417 0.0119 667  -3.510  0.0043 **

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_var_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_var_posthoc.csv')
dt <- read_rds('tables/orc_spr_var_posthoc.rds')
csv <- read_csv('tables/orc_spr_var_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(var ~ dependency * environment * group + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 1 sec
# singular fit

check_convergence(model1) # okay
check_singularity(model1) # warning, but random-effects table looks okay

# set as final model

model <- model1

# save ...

model %>% write_rds('models/src_spr_var_model.rds')
model <- read_rds('models/src_spr_var_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_spr_var_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_var_model.csv')
dt <- read_rds('tables/src_spr_var_model.rds')
csv <- read_csv('tables/src_spr_var_model.csv')

dt %>% rename("Est." = 2) %>% kbl(format = 'markdown')

# |Parameter                                |Est.  |$SE$ |$t$   |$p$    |
# |:----------------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                              |0.06  |0.00 |16.37 |< .001 | ***
# |Dependency[RP]                           |−0.04 |0.00 |−8.67 |< .001 | ***
# |Environment[Long]                        |0.02  |0.01 |3.42  |< .001 | ***
# |Environment[Island]                      |0.03  |0.01 |5.56  |< .001 | ***
# |Group[KLE]                               |0.01  |0.01 |0.65  |0.519  |
# |Group[MLE]                               |0.00  |0.01 |0.12  |0.901  |
# |Dependency[RP]:Environment[Long]         |−0.02 |0.01 |−2.34 |0.019  | *
# |Dependency[RP]:Environment[Island]       |−0.03 |0.01 |−3.21 |0.001  | **
# |Dependency[RP]:Group[KLE]                |0.00  |0.01 |−0.31 |0.756  |
# |Dependency[RP]:Group[MLE]                |0.00  |0.01 |0.24  |0.811  |
# |Environment[Long]:Group[KLE]             |0.01  |0.01 |1.06  |0.291  |
# |Environment[Island]:Group[KLE]           |0.01  |0.01 |1.04  |0.299  |
# |Environment[Long]:Group[MLE]             |0.01  |0.01 |0.56  |0.573  |
# |Environment[Island]:Group[MLE]           |−0.01 |0.01 |−0.60 |0.547  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |−0.01 |0.03 |−0.46 |0.648  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.01  |0.03 |0.43  |0.664  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |−0.04 |0.03 |−1.64 |0.102  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.02  |0.03 |0.83  |0.406  |

# joint tests ...

joint_tests(model)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 192  75.210  <.0001 ***
# environment                    2 192  15.718  <.0001 ***
# group                          2 192   0.240  0.7867
# dependency:environment         2 384   5.511  0.0044 **
# dependency:group               2 192   0.160  0.8522
# environment:group              4 192   0.979  0.4200
# dependency:environment:group   4 384   1.653  0.1604

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df t.ratio p.value
# pronoun - gap short       ENS    -0.0239 0.0139 571  -1.722  0.5501
# pronoun - gap long        ENS    -0.0306 0.0139 571  -2.201  0.2256
# pronoun - gap island      ENS    -0.0691 0.0139 571  -4.972  <.0001 ***
# pronoun - gap short       KLE    -0.0274 0.0135 571  -2.039  0.3183
# pronoun - gap long        KLE    -0.0462 0.0135 571  -3.429  0.0058 **
# pronoun - gap island      KLE    -0.0611 0.0135 571  -4.540  0.0001 ***
# pronoun - gap short       MLE    -0.0141 0.0131 571  -1.081  0.9469
# pronoun - gap long        MLE    -0.0634 0.0131 571  -4.855  <.0001 ***
# pronoun - gap island      MLE    -0.0376 0.0131 571  -2.879  0.0366 *

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% write_rds('tables/src_spr_var_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_var_posthoc.csv')
dt <- read_rds('tables/src_spr_var_posthoc.rds')
csv <- read_csv('tables/src_spr_var_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: accuracy variance proficiency ####
# ---------------------------------------------------------------------------- #
  
# prep dataframe ...

md <- spr_var %>%
  filter(group %in% c('KLE', 'MLE')) %>%
  filter(is.na(proficiency) == FALSE) %>%
  select(study, group, participant, dependency, environment, var, proficiency) %>%
  mutate_at(c('study', 'group', 'participant', 'dependency', 'environment'), factor) %>%
  mutate(var = as.numeric(var)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model2 <- lmer(var ~ dependency * environment * group * proficiency + 
                 (1 + dependency | participant), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 1 sec
# no warnings

check_convergence(model2) # okay
check_singularity(model2) # okay

# set as final model ...

model <- model2

# save ...

model %>% write_rds('models/orc_spr_var_proficiency_model.rds')
model <- read_rds('models/orc_spr_var_proficiency_model.rds')

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3')) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_var_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_var_proficiency_model.csv')
dt <- read_rds('tables/orc_spr_var_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                            |Est.  |$SE$ |df     |$t$   |$p$    |
# |:----------------------------------------------------|:-----|:----|:------|:-----|:------|
# |(Intercept)                                          |0.07  |0.00 |133.00 |16.24 |< .001 | ***
# |Dependency[RP]                                       |−0.05 |0.01 |133.00 |−8.32 |< .001 | ***
# |Environment[Long]                                    |0.00  |0.01 |532.00 |0.21  |0.833  |
# |Environment[Island]                                  |0.00  |0.01 |532.00 |0.46  |0.648  |
# |Group[MLE]                                           |−0.04 |0.01 |133.00 |−4.38 |< .001 | ***
# |Proficiency                                          |−0.03 |0.00 |133.00 |−5.56 |< .001 | ***
# |Dependency[RP]:Environment[Long]                     |0.00  |0.01 |532.00 |0.14  |0.889  |
# |Dependency[RP]:Environment[Island]                   |−0.01 |0.01 |532.00 |−0.83 |0.408  |
# |Dependency[RP]:Group[MLE]                            |0.04  |0.01 |133.00 |3.80  |< .001 | ***
# |Environment[Long]:Group[MLE]                         |0.00  |0.01 |532.00 |−0.38 |0.701  |
# |Environment[Island]:Group[MLE]                       |0.00  |0.01 |532.00 |0.31  |0.756  |
# |Dependency[RP]:Proficiency                           |0.01  |0.01 |133.00 |1.32  |0.188  |
# |Environment[Long]:Proficiency                        |0.00  |0.01 |532.00 |−0.63 |0.531  |
# |Environment[Island]:Proficiency                      |0.00  |0.01 |532.00 |−0.10 |0.920  |
# |Group[MLE]:Proficiency                               |0.01  |0.01 |133.00 |0.96  |0.340  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]               |−0.02 |0.03 |532.00 |−0.84 |0.403  |
# |Dependency[RP]:Envir.[Island]:Group[MLE]             |−0.04 |0.03 |532.00 |−1.48 |0.140  |
# |Dependency[RP]:Envir.[Long]:Proficiency              |0.00  |0.01 |532.00 |−0.36 |0.720  |
# |Dependency[RP]:Envir.[Island]:Proficiency            |0.02  |0.01 |532.00 |1.17  |0.241  |
# |Dependency[RP]:Group[MLE]:Proficiency                |0.00  |0.01 |133.00 |−0.08 |0.938  |
# |Environment[Long]:Group[MLE]:Proficiency             |0.01  |0.01 |532.00 |0.97  |0.331  |
# |Environment[Island]:Group[MLE]:Proficiency           |0.03  |0.01 |532.00 |2.68  |0.008  | **
# |Dependency[RP]:Envir.[Long]:Group[MLE]:Proficiency   |0.04  |0.03 |532.00 |1.58  |0.114  |
# |Dependency[RP]:Envir.[Island]:Group[MLE]:Proficiency |0.02  |0.03 |532.00 |0.69  |0.489  |

# joint tests ...

joint_tests(model)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 133  69.254  <.0001 ***
# environment                                2 532   0.105  0.9006
# group                                      1 133  19.157  <.0001 ***
# proficiency                                1 133  30.920  <.0001 ***
# dependency:environment                     2 532   0.548  0.5785
# dependency:group                           1 133  14.420  0.0002 ***
# dependency:proficiency                     1 133   1.752  0.1879
# environment:group                          2 532   0.242  0.7852
# environment:proficiency                    2 532   0.227  0.7974
# group:proficiency                          1 133   0.916  0.3404
# dependency:environment:group               2 532   1.098  0.3343
# dependency:environment:proficiency         2 532   1.287  0.2771
# dependency:group:proficiency               1 133   0.006  0.9383
# environment:group:proficiency              2 532   3.675  0.0260 *
# dependency:environment:group:proficiency   2 532   1.262  0.2839

pairwise <- model %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment group estimate     SE  df t.ratio p.value
# pronoun - gap short       KLE    0.01395 0.0134 584   1.042  0.8795
# pronoun - gap long        KLE   -0.01116 0.0134 584  -0.833  0.9553
# pronoun - gap island      KLE    0.02017 0.0134 584   1.506  0.5730
# pronoun - gap short       MLE   -0.00651 0.0127 584  -0.511  0.9964
# pronoun - gap long        MLE    0.00932 0.0127 584   0.731  0.9763
# pronoun - gap island      MLE    0.01761 0.0127 584   1.382  0.6661

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_var_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_var_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_spr_var_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(var ~ dependency * environment * group * proficiency + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 1 sec
# singular fit

check_convergence(model1) # okay
check_singularity(model1) # warning, but random-effects table looks okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/src_spr_var_proficiency_model.rds')
model <- read_rds('models/src_spr_var_proficiency_model.rds')

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3')) %>%
  my_table()

dt %>% write_rds('tables/src_spr_var_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_var_proficiency_model.csv')
dt <- read_rds('tables/src_spr_var_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                            |Est.  |$SE$ |df     |$t$   |$p$    |
# |:----------------------------------------------------|:-----|:----|:------|:-----|:------|
# |(Intercept)                                          |0.07  |0.00 |129.42 |14.16 |< .001 | ***
# |Dependency[RP]                                       |−0.04 |0.01 |191.79 |−7.13 |< .001 | ***
# |Environment[Long]                                    |0.02  |0.01 |206.55 |2.93  |0.004  | **
# |Environment[Island]                                  |0.03  |0.01 |407.72 |4.77  |< .001 | ***
# |Group[MLE]                                           |−0.01 |0.01 |129.42 |−1.58 |0.117  |
# |Proficiency                                          |−0.02 |0.00 |129.42 |−4.64 |< .001 | ***
# |Dependency[RP]:Environment[Long]                     |−0.03 |0.01 |515.97 |−2.53 |0.012  | *
# |Dependency[RP]:Environment[Island]                   |−0.02 |0.01 |515.97 |−1.67 |0.096  | .
# |Dependency[RP]:Group[MLE]                            |0.01  |0.01 |191.79 |0.72  |0.470  |
# |Environment[Long]:Group[MLE]                         |−0.01 |0.01 |206.55 |−0.40 |0.693  |
# |Environment[Island]:Group[MLE]                       |−0.02 |0.01 |407.72 |−1.33 |0.186  |
# |Dependency[RP]:Proficiency                           |0.01  |0.01 |191.79 |0.84  |0.399  |
# |Environment[Long]:Proficiency                        |0.00  |0.01 |206.55 |0.52  |0.606  |
# |Environment[Island]:Proficiency                      |0.01  |0.01 |407.72 |0.92  |0.358  |
# |Group[MLE]:Proficiency                               |0.00  |0.01 |129.42 |−0.20 |0.841  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]               |−0.02 |0.03 |515.97 |−0.59 |0.556  |
# |Dependency[RP]:Envir.[Island]:Group[MLE]             |0.01  |0.03 |515.97 |0.44  |0.659  |
# |Dependency[RP]:Envir.[Long]:Proficiency              |0.02  |0.01 |515.97 |1.72  |0.086  | .
# |Dependency[RP]:Envir.[Island]:Proficiency            |0.00  |0.01 |515.97 |−0.07 |0.943  |
# |Dependency[RP]:Group[MLE]:Proficiency                |0.00  |0.01 |191.79 |−0.33 |0.740  |
# |Environment[Long]:Group[MLE]:Proficiency             |−0.02 |0.01 |206.55 |−1.22 |0.224  |
# |Environment[Island]:Group[MLE]:Proficiency           |0.00  |0.01 |407.72 |0.09  |0.925  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]:Proficiency   |−0.02 |0.03 |515.97 |−0.57 |0.567  |
# |Dependency[RP]:Envir.[Island]:Group[MLE]:Proficiency |0.05  |0.03 |515.97 |1.98  |0.048  | *

# joint tests ...

joint_tests(model)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 129  50.864  <.0001 ***
# environment                                2 129  11.630  <.0001 ***
# group                                      1 129   2.486  0.1173
# proficiency                                1 129  21.530  <.0001 ***
# dependency:environment                     2 258   3.313  0.0380 *
# dependency:group                           1 129   0.524  0.4706
# dependency:proficiency                     1 129   0.714  0.3997
# environment:group                          2 129   0.915  0.4030
# environment:proficiency                    2 129   0.427  0.6534
# group:proficiency                          1 129   0.040  0.8410
# dependency:environment:group               2 258   0.534  0.5871
# dependency:environment:proficiency         2 258   2.058  0.1298
# dependency:group:proficiency               1 129   0.110  0.7404
# environment:group:proficiency              2 129   1.041  0.3560
# dependency:environment:group:proficiency   2 258   3.601  0.0287 *

pairwise <- model %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment group estimate     SE  df t.ratio p.value
# pronoun - gap short       KLE    0.00598 0.0144 383   0.414  0.9989
# pronoun - gap long        KLE    0.03629 0.0144 383   2.513  0.0717 .
# pronoun - gap island      KLE   -0.02119 0.0144 383  -1.468  0.6010
# pronoun - gap short       MLE   -0.01042 0.0129 383  -0.806  0.9614
# pronoun - gap long        MLE    0.00475 0.0129 383   0.367  0.9994
# pronoun - gap island      MLE    0.01486 0.0129 383   1.150  0.8210

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% write_rds('tables/src_spr_var_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_var_proficiency_posthoc.csv')
dt <- read_rds('tables/src_spr_var_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: accuracy variance proficiency ####
# ---------------------------------------------------------------------------- #

plot <- md %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = var, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='Variance'),
  coord_cartesian(ylim = c(0, .3)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_spr_var_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_spr_var_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: accuracy filler ####
# ---------------------------------------------------------------------------- #

# check accuracy on grammatical vs. ungrammatical fillers ...

test <- spr_clean %>%
  filter(class == 'filler') %>%
  group_by(study, group, condition) %>%
  summarise(mean = mean(accuracy)) %>%
  ungroup() %>%
  mutate(mean = format(round(mean * 100, 2), nsmall = 2)) %>%
  arrange(condition, study, group) 

kable(test, format = 'markdown')

# summarize ...

plot <- spr_clean %>%
  filter(class == 'filler') %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(accuracy = as.logical(accuracy)) %>%
  ungroup() %>%
  group_by(study, group, group_label, set, condition) %>%
  summarise(mean = mean(accuracy, na.rm=T),
            lcl = my.lcl(accuracy),
            ucl = my.ucl(accuracy)) %>%
  mutate(width = 0.3 * n()) %>%
  ungroup() %>%
  mutate(superset = set) %>%
  mutate(set = case_when(set == 1 & condition == 'grammatical' ~ 'set1',
                         set == 1 & condition == 'ungrammatical' ~ 'set2',
                         set == 2 & condition == 'grammatical' ~ 'set3',
                         set == 2 & condition == 'ungrammatical' ~ 'set4',
                         set == 3 ~ 'set5',
                         set == 4 ~ 'set6'))

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = superset, y = mean * 100, col = condition, shape = condition, width = width),
  geom_line(aes(linetype = condition), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl * 100, ymax = ucl * 100), 
                linetype = 1, position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(50, 100)),
  scale_x_discrete(name="Filler Type", 
                   limits = c('1', '2', '3', '4'),
                   labels = c('F1', 'F2', 'F3', 'F4')),
  scale_y_continuous(name="Accuracy (%)", breaks=c(0, 25, 50, 75, 100)),
  scale_color_manual(name = 'Grammaticality', 
                     values = c('#40B0A6', 'coral'),
                     labels=c('Grammatical', 'Ungrammatical')),
  scale_linetype_manual(name = 'Grammaticality', 
                     values=c('solid', '41'), 
                     labels=c('Grammatical', 'Ungrammatical')),
  theme(text = element_text(size = 12), 
        axis.title.y = element_text(margin=margin(r=-3)),
        axis.text.x = element_markdown(),
        strip.text = element_markdown(),
        legend.margin=margin(0, 0, 0, -5)),
  facet_wrap(~group_label)
)

p1 + s
ggsave("plots/orc_spr_accuracy_filler.png", width=6.5, height=2.25, dpi=600)

p2 + s
ggsave("plots/src_spr_accuracy_filler.png", width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot %>%
  mutate(mean = mean * 100,
         lcl = lcl * 100,
         ucl = ucl * 100) %>%
  select(-group_label, -width, -superset) %>%
  rename( 'Type' = 'set',
         'Accuracy (\\%)' = 'mean') %>%
  my_table() 

dt <- dt %>%
  mutate(Type = case_when(Type == "F1" ~ "F1G",
                          Type == "F2" ~ "F1U",
                          Type == "F3" ~ "F2G",
                          Type == "F4" ~ "F2U",
                          Type == "F5" ~ "F3G",
                          Type == "F6" ~ "F4G")) %>%
  rename("Grammaticality" = "Condition") 

dt %>% write_rds('tables/spr_accuracy_filler.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/spr_accuracy_filler.csv')
dt <- read_rds('tables/spr_accuracy_filler.rds')

# ---------------------------------------------------------------------------- #
# model: accuracy ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_clean %>%
  filter(class == 'critical') %>%
  group_by(study, group, participant) %>%
  mutate(acc_rate = mean(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(study, group, participant, item, dependency, environment, acc_rate, accuracy) %>%
  summarise() %>%
  ungroup() %>%
  mutate(dependency = fct_drop(dependency),
         environment = fct_drop(environment),
         accuracy = as.logical(accuracy)) %>%
  mutate_at(c('study', 'group', 'participant', 'item'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

# apply deviation coding ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

# view contrasts ...

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4b <- glmer(accuracy ~ dependency * environment * group + 
                  (1 + dependency + environment | participant) + 
                  (1 + dependency + environment + group | item), 
                data = filter(md, study == 'ORC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model4b)
toc()
beep()

# ... 11350.9 sec elapsed

performance::check_convergence(model4b) # convergence: okay
performance::check_singularity(model4b) # singularity: warning

# save ...

model4b %>% write_rds('models/orc_spr_accuracy_model4b.rds')
model4b <- read_rds('models/orc_spr_accuracy_model4b.rds')

# summary ...

summary(model4b)

# Parameter                        Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                      3.087614   0.142592  21.654  < 2e-16 ***
# dependency2                      1.012369   0.173932   5.820 5.87e-09 ***
# environment2                    -0.067325   0.220390  -0.305 0.759999    
# environment3                    -0.152502   0.219653  -0.694 0.487504    
# group2                          -1.721434   0.267425  -6.437 1.22e-10 ***
# group3                          -0.929063   0.274962  -3.379 0.000728 ***
# dependency2:environment2         0.015119   0.277042   0.055 0.956480    
# dependency2:environment3         0.080574   0.279485   0.288 0.773122    
# dependency2:group2               0.653370   0.276534   2.363 0.018142 *  
# dependency2:group3               0.029023   0.287447   0.101 0.919576    
# environment2:group2             -0.243852   0.349879  -0.697 0.485826    
# environment3:group2              0.002839   0.350983   0.008 0.993547    
# environment2:group3             -0.263678   0.363899  -0.725 0.468704    
# environment3:group3             -0.163385   0.367145  -0.445 0.656310    
# dependency2:environment2:group2 -0.172552   0.646544  -0.267 0.789559    
# dependency2:environment3:group2 -0.187235   0.639293  -0.293 0.769616    
# dependency2:environment2:group3 -0.237029   0.683124  -0.347 0.728607    
# dependency2:environment3:group3  0.280383   0.675291   0.415 0.677993 

# table ...

dt <- summary(model4b)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/orc_spr_accuracy_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_accuracy_model.csv')
dt <- read_rds('tables/orc_spr_accuracy_model.rds')

# joint tests ...

joint_tests(model4b)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  33.878  <.0001
# environment                    2 Inf   0.251  0.7784
# group                          2 Inf  20.721  <.0001
# dependency:environment         2 Inf   0.049  0.9522
# dependency:group               2 Inf   4.280  0.0138
# environment:group              4 Inf   0.256  0.9061
# dependency:environment:group   4 Inf   0.311  0.8711

# post-hoc tests ...

pairwise <- model4b %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate    SE  df z.ratio p.value
# pronoun - gap short       ENS      0.718 0.411 Inf   1.745  0.5139
# pronoun - gap long        ENS      0.869 0.416 Inf   2.092  0.2735
# pronoun - gap island      ENS      0.767 0.399 Inf   1.922  0.3822
# pronoun - gap short       KLE      1.491 0.300 Inf   4.965  <.0001
# pronoun - gap long        KLE      1.470 0.281 Inf   5.231  <.0001
# pronoun - gap island      KLE      1.353 0.286 Inf   4.726  <.0001
# pronoun - gap short       MLE      0.732 0.360 Inf   2.034  0.3080
# pronoun - gap long        MLE      0.647 0.340 Inf   1.906  0.3937
# pronoun - gap island      MLE      1.062 0.339 Inf   3.133  0.0153

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_accuracy_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_accuracy_posthoc.csv')
dt <- read_rds('tables/orc_spr_accuracy_posthoc.rds')

# plot ...

plot(emmeans(model4b, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4b <- glmer(accuracy ~ dependency * environment * group + 
                  (1 + dependency + environment | participant) + 
                  (1 + dependency + environment + group | item), 
                data = filter(md, study == 'SRC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model4b)
toc()
beep()

# ... 3435.53 sec elapsed (57 min)

performance::check_convergence(model4b) # convergence: warning (but not lme4)
performance::check_singularity(model4b) # singularity: warning

# save ...

model4b %>% write_rds('models/src_spr_accuracy_model4b.rds')
model4b <- read_rds('models/src_spr_accuracy_model4b.rds')

# model summary...

summary(model4b)

# Parameter                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                      2.86555    0.14375  19.934  < 2e-16 ***
# dependency2                      1.13288    0.16543   6.848 7.49e-12 ***
# environment2                    -0.63791    0.23399  -2.726  0.00641 ** 
# environment3                    -0.92247    0.20972  -4.399 1.09e-05 ***
# group2                          -0.32159    0.26018  -1.236  0.21646    
# group3                          -0.16969    0.26711  -0.635  0.52525    
# dependency2:environment2         0.37387    0.28752   1.300  0.19348    
# dependency2:environment3         0.51112    0.27704   1.845  0.06504 .  
# dependency2:group2               0.17276    0.28430   0.608  0.54341    
# dependency2:group3               0.09463    0.28542   0.332  0.74024    
# environment2:group2             -0.56249    0.36004  -1.562  0.11822    
# environment3:group2             -0.34668    0.34159  -1.015  0.31015    
# environment2:group3             -0.33838    0.36459  -0.928  0.35335    
# environment3:group3             -0.19883    0.33975  -0.585  0.55840    
# dependency2:environment2:group2  0.43253    0.65319   0.662  0.50786    
# dependency2:environment3:group2 -0.30589    0.64336  -0.475  0.63447    
# dependency2:environment2:group3  1.14944    0.66000   1.742  0.08158 .  
# dependency2:environment3:group3 -0.46426    0.63736  -0.728  0.46636 

# table ...

dt <- summary(model4b)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_spr_accuracy_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_accuracy_model.csv')
dt <- read_rds('tables/src_spr_accuracy_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                     | Estimate|   SE|     z|     p|
# |:---------------------------------------------|--------:|----:|-----:|-----:|
# |(Intercept)                                   |     2.87| 0.14| 19.93| 0.000| ***
# |Dependency[RP]                                |     1.13| 0.17|  6.85| 0.000| ***
# |Environment[Long]                             |    -0.64| 0.23| -2.73| 0.006| **
# |Environment[Island]                           |    -0.92| 0.21| -4.40| 0.000| ***
# |Group[KLE]                                    |    -0.32| 0.26| -1.24| 0.216|
# |Group[MLE]                                    |    -0.17| 0.27| -0.64| 0.525|
# |Dependency[RP]:Environment[Long]              |     0.37| 0.29|  1.30| 0.193|
# |Dependency[RP]:Environment[Island]            |     0.51| 0.28|  1.84| 0.065| .
# |Dependency[RP]:Group[KLE]                     |     0.17| 0.28|  0.61| 0.543|
# |Dependency[RP]:Group[MLE]                     |     0.09| 0.29|  0.33| 0.740|
# |Environment[Long]:Group[KLE]                  |    -0.56| 0.36| -1.56| 0.118|
# |Environment[Island]:Group[KLE]                |    -0.35| 0.34| -1.01| 0.310|
# |Environment[Long]:Group[MLE]                  |    -0.34| 0.36| -0.93| 0.353|
# |Environment[Island]:Group[MLE]                |    -0.20| 0.34| -0.59| 0.558|
# |Dependency[RP]:Environment[Long]:Group[KLE]   |     0.43| 0.65|  0.66| 0.508|
# |Dependency[RP]:Environment[Island]:Group[KLE] |    -0.31| 0.64| -0.48| 0.634|
# |Dependency[RP]:Environment[Long]:Group[MLE]   |     1.15| 0.66|  1.74| 0.082| .
# |Dependency[RP]:Environment[Island]:Group[MLE] |    -0.46| 0.64| -0.73| 0.466|


# joint tests ...

joint_tests(model4b)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  46.894  <.0001
# environment                    2 Inf   9.768  0.0001
# group                          2 Inf   0.767  0.4646
# dependency:environment         2 Inf   1.733  0.1767
# dependency:group               2 Inf   0.185  0.8314
# environment:group              4 Inf   0.612  0.6539
# dependency:environment:group   4 Inf   1.907  0.1062

# post-hoc tests ...

pairwise <- model4b %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate    SE  df z.ratio p.value
# pronoun - gap short       ENS      0.839 0.399 Inf   2.102  0.2705
# pronoun - gap long        ENS      0.686 0.358 Inf   1.916  0.3901
# pronoun - gap island      ENS      1.607 0.344 Inf   4.667  <.0001 ***
# pronoun - gap short       KLE      0.970 0.395 Inf   2.458  0.1165
# pronoun - gap long        KLE      1.249 0.313 Inf   3.985  0.0006 ***
# pronoun - gap island      KLE      1.431 0.304 Inf   4.706  <.0001 ***
# pronoun - gap short       MLE      0.705 0.371 Inf   1.900  0.4013
# pronoun - gap long        MLE      1.701 0.336 Inf   5.070  <.0001 ***
# pronoun - gap island      MLE      1.009 0.296 Inf   3.410  0.0058 **

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_accuracy_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_accuracy_posthoc.csv')
dt <- read_rds('tables/src_spr_accuracy_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|    z|     p|
# |:-----|:-----------|:--------|--------:|----:|----:|-----:|
# |ENS   |Short       |RP - Gap |     0.84| 0.40| 2.10| 0.270|
# |ENS   |Long        |RP - Gap |     0.69| 0.36| 1.92| 0.390|
# |ENS   |Island      |RP - Gap |     1.61| 0.34| 4.67| 0.000| ***
# |KLE   |Short       |RP - Gap |     0.97| 0.39| 2.46| 0.117|
# |KLE   |Long        |RP - Gap |     1.25| 0.31| 3.99| 0.001| ***
# |KLE   |Island      |RP - Gap |     1.43| 0.30| 4.71| 0.000| ***
# |MLE   |Short       |RP - Gap |     0.71| 0.37| 1.90| 0.401|
# |MLE   |Long        |RP - Gap |     1.70| 0.34| 5.07| 0.000| ***
# |MLE   |Island      |RP - Gap |     1.01| 0.30| 3.41| 0.006| **


# plot ...

plot(emmeans(model4b, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: accuracy proficiency ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_clean %>%
  filter(class == 'critical') %>%
  filter(group %in% c('KLE', 'MLE')) %>%
  filter(is.na(proficiency) == FALSE) %>%
  group_by(study, group, participant) %>%
  mutate(acc_rate = mean(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(study, group, participant, item, dependency, environment, acc_rate, accuracy, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  mutate(dependency = fct_drop(dependency),
         environment = fct_drop(environment),
         accuracy = as.logical(accuracy)) %>%
  mutate_at(c('study', 'group', 'participant', 'item'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

# check participants ...

check <- md %>%
  group_by(study, group) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant, proficiency, dependency, environment) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = accuracy * 100, col = dependency, fill = dependency),
  # geom_hline(yintercept = 50),
  # geom_vline(xintercept = 0),
  geom_smooth(aes(linetype = dependency), 
              method=lm, 
              alpha = .75,
              linewidth = .5),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name="Accuracy (%)",
                     breaks = c(0, 25, 50, 75, 100)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                    values=c('41', 'solid'), 
                    labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        axis.title.x = element_markdown(),
        strip.text.y = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(50, 100))
ggsave('plots/orc_spr_accuracy_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(50, 100))
ggsave('plots/src_spr_accuracy_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4b <- glmer(accuracy ~ dependency * environment * group * proficiency + 
                  (1 + dependency + environment | participant) + 
                  (1 + dependency + environment + group | item), 
                data = filter(md, study == 'ORC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model4b)
toc()
beep()

# ... 1703.73 sec elapsed

performance::check_convergence(model4b) # warning (but not lme4)
performance::check_singularity(model4b) # warning

model4b %>% write_rds('models/orc_spr_accuracy_proficiency_model4b.rds')
model4b <- read_rds('models/orc_spr_accuracy_proficiency_model4b.rds')

# model summary...

summary(model4b)

# Parameter                                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  2.712398   0.160838  16.864  < 2e-16 ***
# dependency2                                  1.271595   0.190072   6.690 2.23e-11 ***
# environment2                                -0.149656   0.237753  -0.629  0.52905    
# environment3                                -0.240757   0.236561  -1.018  0.30880    
# group2                                       0.827309   0.258140   3.205  0.00135 ** 
# proficiency                                  0.716151   0.123843   5.783 7.35e-09 ***
# dependency2:environment2                    -0.001501   0.328216  -0.005  0.99635    
# dependency2:environment3                     0.095705   0.327447   0.292  0.77008    
# dependency2:group2                          -0.633865   0.265882  -2.384  0.01713 *  
# environment2:group2                         -0.166234   0.348161  -0.477  0.63303    
# environment3:group2                         -0.423863   0.340236  -1.246  0.21284    
# dependency2:proficiency                     -0.036453   0.118557  -0.307  0.75848    
# environment2:proficiency                     0.114826   0.155095   0.740  0.45908    
# environment3:proficiency                     0.054664   0.152726   0.358  0.72040    
# group2:proficiency                           0.128742   0.245567   0.524  0.60010    
# dependency2:environment2:group2             -0.027132   0.613575  -0.044  0.96473    
# dependency2:environment3:group2              0.374595   0.610212   0.614  0.53930    
# dependency2:environment2:proficiency         0.206281   0.278033   0.742  0.45813    
# dependency2:environment3:proficiency         0.028588   0.278576   0.103  0.91826    
# dependency2:group2:proficiency              -0.209633   0.233314  -0.898  0.36892    
# environment2:group2:proficiency             -0.507154   0.300141  -1.690  0.09108 .  
# environment3:group2:proficiency             -0.624229   0.301802  -2.068  0.03861 *  
# dependency2:environment2:group2:proficiency  0.205701   0.553813   0.371  0.71032    
# dependency2:environment3:group2:proficiency -0.080931   0.556364  -0.145  0.88434 


# table ...

dt <- summary(model4b)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'Group[MLE]')) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_accuracy_proficiency_model.rds')

# joint tests ...

joint_tests(model4b)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  44.757  <.0001
# environment                                2 Inf   0.520  0.5946
# group                                      1 Inf  10.271  0.0014
# proficiency                                1 Inf  33.440  <.0001
# dependency:environment                     2 Inf   0.062  0.9398
# dependency:group                           1 Inf   5.683  0.0171
# dependency:proficiency                     1 Inf   0.095  0.7585
# environment:group                          2 Inf   0.800  0.4492
# environment:proficiency                    2 Inf   0.275  0.7599
# group:proficiency                          1 Inf   0.275  0.6001
# dependency:environment:group               2 Inf   0.279  0.7566
# dependency:environment:proficiency         2 Inf   0.332  0.7178
# dependency:group:proficiency               1 Inf   0.807  0.3689
# environment:group:proficiency              2 Inf   2.391  0.0916
# dependency:environment:group:proficiency   2 Inf   0.151  0.8602

# post-hoc tests ...

pairwise <- model4b %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate    SE  df z.ratio p.value
# pronoun - gap short       KLE     0.0109 0.273 Inf   0.040  1.0000
# pronoun - gap long        KLE     0.1143 0.261 Inf   0.437  0.9985
# pronoun - gap island      KLE     0.0799 0.266 Inf   0.301  0.9998
# pronoun - gap short       MLE    -0.2404 0.306 Inf  -0.787  0.9662
# pronoun - gap long        MLE     0.0688 0.285 Inf   0.241  1.0000
# pronoun - gap island      MLE    -0.2522 0.285 Inf  -0.885  0.9409

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

# dt <- pairwise %>%
#   rename('Group' = 'group', 'Environment' = 'environment', 
#          'Contrast' = 'contrast', 'z' = 'z.ratio', 'Estimate' = 'estimate',
#          'p' = 'p.value') %>%
#   select(Group, Environment, Contrast, everything()) %>%
#   mutate_at(c('Environment', 'Contrast'), as.character) %>%
#   mutate_if(is.character, str_replace_all, 
#             c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
#               'gap' = 'Gap', 'pronoun' = 'RP')) %>%
#   select(-df) %>%
#   mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
#   arrange(Group, Environment) %>%
#   mutate(across(4:6, round, 2), across(7, round, 3))

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast | Estimate|   SE|     z|     p|
# |:-----|:-----------|:--------|--------:|----:|-----:|-----:|
# |KLE   |Short       |RP - Gap |     0.01| 0.27|  0.04| 1.000|
# |KLE   |Long        |RP - Gap |     0.11| 0.26|  0.44| 0.999|
# |KLE   |Island      |RP - Gap |     0.08| 0.27|  0.30| 1.000|
# |MLE   |Short       |RP - Gap |    -0.24| 0.31| -0.79| 0.966|
# |MLE   |Long        |RP - Gap |     0.07| 0.29|  0.24| 1.000|
# |MLE   |Island      |RP - Gap |    -0.25| 0.28| -0.89| 0.941|

dt %>% write_rds('tables/orc_spr_accuracy_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_accuracy_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_spr_accuracy_proficiency_posthoc.rds')

# plot ...

emmip(model4b, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model5b <- glmer(accuracy ~ dependency * environment * group * proficiency + 
                   (1 + dependency | participant) + 
                   (1 + dependency + group | item), 
                 data = filter(md, study == 'SRC'), family = binomial, 
                 control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model5b)
toc()
beep()

# ... 121.04 sec elapsed

performance::check_convergence(model5b) # okay
performance::check_singularity(model5b) # warning

# save ...

model5b %>% write_rds('models/src_spr_accuracy_proficiency_model5b.rds')
model5b <- read_rds('models/src_spr_accuracy_proficiency_model5b.rds')

# model summary...

summary(model5b)

# Parameter                                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  2.74598    0.15101  18.185  < 2e-16 ***
# dependency2                                  1.12842    0.18269   6.177 6.54e-10 ***
# environment2                                -0.66542    0.17745  -3.750 0.000177 ***
# environment3                                -0.93862    0.16224  -5.785 7.24e-09 ***
# group2                                       0.41583    0.24202   1.718 0.085772 .  
# proficiency                                  0.50630    0.11562   4.379 1.19e-05 ***
# dependency2:environment2                     0.88402    0.35512   2.489 0.012798 *  
# dependency2:environment3                     0.44172    0.32584   1.356 0.175217    
# dependency2:group2                           0.17938    0.27980   0.641 0.521452    
# environment2:group2                          0.36067    0.35489   1.016 0.309491    
# environment3:group2                         -0.00612    0.32500  -0.019 0.984976    
# dependency2:proficiency                      0.19789    0.12724   1.555 0.119877    
# environment2:proficiency                    -0.20636    0.15947  -1.294 0.195650    
# environment3:proficiency                    -0.31737    0.14840  -2.139 0.032466 *  
# group2:proficiency                           0.03470    0.23059   0.150 0.880393    
# dependency2:environment2:group2              0.75465    0.70967   1.063 0.287609    
# dependency2:environment3:group2             -0.41992    0.65096  -0.645 0.518881    
# dependency2:environment2:proficiency        -0.23174    0.31737  -0.730 0.465267    
# dependency2:environment3:proficiency        -0.06579    0.29547  -0.223 0.823802    
# dependency2:group2:proficiency               0.32729    0.25286   1.294 0.195535    
# environment2:group2:proficiency              0.46872    0.32022   1.464 0.143262    
# environment3:group2:proficiency             -0.04352    0.29786  -0.146 0.883829    
# dependency2:environment2:group2:proficiency  0.80423    0.63601   1.264 0.206055    
# dependency2:environment3:group2:proficiency -0.29046    0.59244  -0.490 0.623943 

# table ...

dt <- summary(model5b)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_spr_accuracy_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_accuracy_proficiency_model.csv')
dt <- read_rds('tables/src_spr_accuracy_proficiency_model.rds')

# joint tests ...

joint_tests(model5b)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  38.153  <.0001
# environment                                2 Inf  16.775  <.0001
# group                                      1 Inf   2.952  0.0858
# proficiency                                1 Inf  19.176  <.0001
# dependency:environment                     2 Inf   3.120  0.0442
# dependency:group                           1 Inf   0.411  0.5215
# dependency:proficiency                     1 Inf   2.419  0.1199
# environment:group                          2 Inf   0.839  0.4321
# environment:proficiency                    2 Inf   2.288  0.1015
# group:proficiency                          1 Inf   0.023  0.8804
# dependency:environment:group               2 Inf   1.892  0.1508
# dependency:environment:proficiency         2 Inf   0.299  0.7414
# dependency:group:proficiency               1 Inf   1.675  0.1955
# environment:group:proficiency              2 Inf   1.849  0.1574
# dependency:environment:group:proficiency   2 Inf   1.956  0.1415

# post-hoc tests ...

pairwise <- model5b %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate    SE  df z.ratio p.value
# pronoun - gap short       KLE      0.219 0.330 Inf   0.663  0.9857
# pronoun - gap long        KLE     -0.415 0.278 Inf  -1.492  0.5828
# pronoun - gap island      KLE      0.298 0.262 Inf   1.141  0.8274
# pronoun - gap short       MLE      0.375 0.344 Inf   1.090  0.8552
# pronoun - gap long        MLE      0.545 0.336 Inf   1.625  0.4830
# pronoun - gap island      MLE      0.164 0.250 Inf   0.656  0.9864

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(group, environment) %>%
  my_table()

dt %>% kable(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$   |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:-----|
# |KLE   |Short       |RP − Gap |0.22                     |0.33 |0.66  |0.986 |
# |KLE   |Long        |RP − Gap |−0.41                    |0.28 |−1.49 |0.583 |
# |KLE   |Island      |RP − Gap |0.30                     |0.26 |1.14  |0.827 |
# |MLE   |Short       |RP − Gap |0.38                     |0.34 |1.09  |0.855 |
# |MLE   |Long        |RP − Gap |0.55                     |0.34 |1.62  |0.483 |
# |MLE   |Island      |RP − Gap |0.16                     |0.25 |0.66  |0.986 |

dt %>% write_rds('tables/src_spr_accuracy_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_accuracy_proficiency_posthoc.csv')
dt <- read_rds('tables/src_spr_accuracy_proficiency_posthoc.rds')

# plot ...

emmip(model5b, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: question RT ####
# ---------------------------------------------------------------------------- #

# summarize to comprehension question data ...

temp1 <- spr_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant, proficiency, item, dependency, environment, question_rt) %>%
  summarise() %>%
  ungroup()

# remove extreme values ...

temp2 <- temp1 %>%
  filter(question_rt <= 30000, question_rt >= 100)

# calculate percent excluded ...

check <- temp1 %>%
  mutate(keep = case_when(question_rt <= 30000 & question_rt >= 100 ~ TRUE,
                          TRUE ~ FALSE))

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check1 %>% kbl(format = 'markdown')

# |study |      prop|
# |:-----|---------:|
# |ORC   | 0.8223201|
# |SRC   | 1.4017094|

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 %>% kbl(format = 'markdown')

# |study |group |      prop|
# |:-----|:-----|---------:|
# |ORC   |ENS   | 0.4074074|
# |ORC   |KLE   | 0.8457711|
# |ORC   |MLE   | 1.3333333|
# |SRC   |ENS   | 1.0928962|
# |SRC   |KLE   | 0.6666667|
# |SRC   |MLE   | 2.3671498|

# remove outliers more than 2 SDs from the mean ...

spr_qrt <- temp2 %>%
  group_by(study, group, dependency, environment) %>%
  mutate(upper = mean(question_rt) + (2 * sd(question_rt)),
         lower = mean(question_rt) - (2 * sd(question_rt))) %>%
  ungroup() %>%
  filter(question_rt <= upper & question_rt >= lower) %>%
  select(-upper, -lower)

# calculate percent excluded ...

check <- temp2 %>%
  group_by(study, group, dependency, environment) %>%
  mutate(upper = mean(question_rt) + (2 * sd(question_rt)),
         lower = mean(question_rt) - (2 * sd(question_rt))) %>%
  ungroup() %>%
  mutate(keep = case_when(question_rt <= upper & question_rt >= lower ~ TRUE,
                          TRUE ~ FALSE)) %>%
  select(study, group, participant, proficiency, dependency, environment, question_rt, upper, lower, keep)

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check1 %>% kbl(format = "markdown")

# |study |     prop|
# |:-----|--------:|
# |ORC   | 5.137696|
# |SRC   | 4.906380|

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 %>% kbl(format = "markdown")

# |study |group |     prop|
# |:-----|:-----|--------:|
# |ORC   |ENS   | 5.355151|
# |ORC   |KLE   | 4.917210|
# |ORC   |MLE   | 5.067568|
# |SRC   |ENS   | 4.917127|
# |SRC   |KLE   | 4.543108|
# |SRC   |MLE   | 5.244928|

# plot (overall) ...

plot <- spr_qrt %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean = mean(question_rt, na.rm=T),
            lcl = my.lcl(question_rt),
            ucl = my.ucl(question_rt),
            mean2 = mean(log(question_rt), na.rm=T),
            lcl2 = my.lcl(log(question_rt)),
            ucl2 = my.ucl(log(question_rt))) %>%
  ungroup()

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

d1 <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE)
)

d2 <- list(
  aes(x=environment, y=mean2, group=dependency, col=dependency, shape=dependency),
  geom_errorbar(aes(ymin=lcl2, ymax=ucl2), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE)
)

s <- list(
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Environment", 
                   limits = c("short", "long", "island"), 
                   labels = c("Short", "Long", "Island")),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    strip.text = element_markdown(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin=margin(0, 0, 0, -5)),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

p1 + d1 + s + 
  coord_cartesian(ylim = c(4000, 7000)) +
  scale_y_continuous(name="RT (ms)", breaks = c(4000, 5000, 6000, 7000))
ggsave("plots/orc_spr_question_rt_overall.png", width=6.5, height=2, dpi=600)

p2 + d1 + s + 
  coord_cartesian(ylim = c(4000, 7000)) +
  scale_y_continuous(name="RT (ms)", breaks = c(4000, 5000, 6000, 7000))
ggsave("plots/src_spr_question_rt_overall.png", width=6.5, height=2, dpi=600)

p1 + d2 + s + 
  coord_cartesian(ylim = c(8.25, 8.75)) +
  scale_y_continuous(name="Log RT", breaks = c(8.25, 8.5, 8.75))
ggsave("plots/orc_spr_question_logrt_overall.png", width=6.5, height=2, dpi=600)

p2 + d2 + s + 
  coord_cartesian(ylim = c(8.25, 8.75)) +
  scale_y_continuous(name="Log RT", breaks = c(8.25, 8.5, 8.75))
ggsave("plots/src_spr_question_logrt_overall.png", width=6.5, height=2, dpi=600)

# raw rt and log rt tables ...

dt <- plot %>%
  select(study, group, dependency, environment, mean, lcl, ucl) %>%
  my_table() %>%
  csv_table()

dt %>% write_csv("tables/spr_question_rt.csv")

dt <- plot %>%
  select(study, group, dependency, environment, mean2, lcl2, ucl2) %>%
  rename(mean = mean2, lcl = lcl2, ucl = ucl2) %>%
  my_table() %>%
  csv_table()

dt %>% write_csv("tables/spr_question_logrt.csv")

# ---------------------------------------------------------------------------- #
# model: question RT ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_qrt %>%
  select(study, group, participant, item, dependency, environment, question_rt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model2 <- lmer(log(question_rt) ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 409.12 sec elapsed

check_convergence(model2) # okay
check_singularity(model2) # warning (but random effects table looks okay)

# save ...

model2 %>% write_rds('models/orc_spr_question_logrt_model2.rds')

# summary ...

summary(model2)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                      8.449e+00  1.678e-02  2.136e+02 503.487  < 2e-16 ***
# dependency2                     -7.587e-02  9.773e-03  3.351e+01  -7.764 5.45e-09 ***
# environment2                     2.412e-02  1.074e-02  4.714e+01   2.245   0.0295 *  
# environment3                     2.325e-02  1.122e-02  3.568e+01   2.072   0.0456 *  
# group2                          -3.028e-02  3.728e-02  2.174e+02  -0.812   0.4176    
# group3                           9.617e-02  3.772e-02  2.236e+02   2.550   0.0115 *  
# dependency2:environment2        -3.057e-02  1.788e-02  5.675e+03  -1.710   0.0874 .  
# dependency2:environment3        -2.065e-03  1.787e-02  5.602e+03  -0.116   0.9080    
# dependency2:group2              -4.194e-02  1.866e-02  2.695e+02  -2.247   0.0254 *  
# dependency2:group3              -4.134e-02  1.848e-02  2.710e+02  -2.237   0.0261 *  
# environment2:group2              2.043e-02  2.334e-02  2.753e+02   0.875   0.3822    
# environment3:group2              9.212e-03  2.261e-02  3.935e+02   0.407   0.6839    
# environment2:group3              6.076e-02  2.317e-02  2.792e+02   2.623   0.0092 ** 
# environment3:group3              3.828e-02  2.251e-02  3.985e+02   1.701   0.0898 .  
# dependency2:environment2:group2  5.527e-02  4.283e-02  5.656e+03   1.290   0.1970    
# dependency2:environment3:group2  4.868e-02  4.286e-02  5.671e+03   1.136   0.2561    
# dependency2:environment2:group3  6.178e-03  4.238e-02  5.690e+03   0.146   0.8841    
# dependency2:environment3:group3  3.043e-03  4.249e-02  5.700e+03   0.072   0.9429

# table ...

dt <- summary(model2)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_spr_question_logrt_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_question_logrt_model.csv')
dt <- read_rds('tables/orc_spr_question_logrt_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                |Estimate|$SE$ |$t$    |$p$    |
# |:----------------------------------------|:-------|:----|:------|:------|
# |(Intercept)                              |8.45    |0.02 |503.49 |< .001 | ***
# |Dependency[RP]                           |−0.08   |0.01 |−7.76  |< .001 | ***
# |Environment[Long]                        |0.02    |0.01 |2.25   |0.029  | *
# |Environment[Island]                      |0.02    |0.01 |2.07   |0.046  | *
# |Group[KLE]                               |−0.03   |0.04 |−0.81  |0.418  |
# |Group[MLE]                               |0.10    |0.04 |2.55   |0.011  | *
# |Dependency[RP]:Environment[Long]         |−0.03   |0.02 |−1.71  |0.087  | .
# |Dependency[RP]:Environment[Island]       |0.00    |0.02 |−0.12  |0.908  |
# |Dependency[RP]:Group[KLE]                |−0.04   |0.02 |−2.25  |0.025  | *
# |Dependency[RP]:Group[MLE]                |−0.04   |0.02 |−2.24  |0.026  | *
# |Environment[Long]:Group[KLE]             |0.02    |0.02 |0.88   |0.382  |
# |Environment[Island]:Group[KLE]           |0.01    |0.02 |0.41   |0.684  |
# |Environment[Long]:Group[MLE]             |0.06    |0.02 |2.62   |0.009  | **
# |Environment[Island]:Group[MLE]           |0.04    |0.02 |1.70   |0.090  | .
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |0.06    |0.04 |1.29   |0.197  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.05    |0.04 |1.14   |0.256  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |0.01    |0.04 |0.15   |0.884  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.00    |0.04 |0.07   |0.943  |

# joint tests ...

joint_tests(model2)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  60.279  <.0001 ***
# environment                    2 Inf   3.257  0.0385 *
# group                          2 Inf   5.441  0.0043 **
# dependency:environment         2 Inf   1.828  0.1608
# dependency:group               2 Inf   3.511  0.0299 *
# environment:group              4 Inf   1.887  0.1096
# dependency:environment:group   4 Inf   0.566  0.6874

# post-hoc tests ...

pairwise <- model2 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS    -0.0247 0.0211 Inf  -1.166  0.9158
# pronoun - gap long        ENS    -0.0757 0.0211 Inf  -3.592  0.0029 **
# pronoun - gap island      ENS    -0.0440 0.0211 Inf  -2.081  0.2877
# pronoun - gap short       KLE    -0.1012 0.0243 Inf  -4.173  0.0003 ***
# pronoun - gap long        KLE    -0.0970 0.0242 Inf  -4.004  0.0006 ***
# pronoun - gap island      KLE    -0.0719 0.0242 Inf  -2.969  0.0265 *
# pronoun - gap short       MLE    -0.0691 0.0238 Inf  -2.908  0.0321 *
# pronoun - gap long        MLE    -0.1140 0.0238 Inf  -4.789  <.0001 ***
# pronoun - gap island      MLE    -0.0853 0.0239 Inf  -3.564  0.0033 **

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_question_logrt_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_question_logrt_posthoc.csv')
dt <- read_rds('tables/orc_spr_question_logrt_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:------|
# |ENS   |Short       |RP − Gap |−0.02                    |0.02 |−1.17 |0.916  |
# |ENS   |Long        |RP − Gap |−0.08                    |0.02 |−3.59 |0.003  | *
# |ENS   |Island      |RP − Gap |−0.04                    |0.02 |−2.08 |0.288  |
# |KLE   |Short       |RP − Gap |−0.10                    |0.02 |−4.17 |< .001 | ***
# |KLE   |Long        |RP − Gap |−0.10                    |0.02 |−4.00 |< .001 | ***
# |KLE   |Island      |RP − Gap |−0.07                    |0.02 |−2.97 |0.026  | *
# |MLE   |Short       |RP − Gap |−0.07                    |0.02 |−2.91 |0.032  | *
# |MLE   |Long        |RP − Gap |−0.11                    |0.02 |−4.79 |< .001 | ***
# |MLE   |Island      |RP − Gap |−0.09                    |0.02 |−3.56 |0.003  | **

# ---------------------------------------------------------------------------- #
# ... src model 1 (final) ####
# ---------------------------------------------------------------------------- #

tic()
model1 <- lmer(log(question_rt) ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 485.52 sec elapsed

check_convergence(model1) # okay
check_singularity(model1) # warning (but random effects table looks okay)

# save ...

model1 %>% write_rds('models/src_spr_question_logrt_model.rds')
model1 <- read_rds('models/src_spr_question_logrt_model.rds')

# summary ...

summary(model1)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                      8.466e+00  1.847e-02  1.915e+02 458.384  < 2e-16 ***
# dependency2                     -7.266e-02  9.164e-03  1.301e+02  -7.929  8.8e-13 ***
# environment2                     1.088e-02  1.036e-02  2.057e+02   1.050  0.29478    
# environment3                     1.110e-01  1.398e-02  4.518e+01   7.940  4.1e-10 ***
# group2                          -1.040e-01  4.365e-02  1.834e+02  -2.382  0.01822 *  
# group3                           1.229e-01  4.363e-02  1.863e+02   2.816  0.00539 ** 
# dependency2:environment2        -2.194e-02  1.906e-02  4.841e+03  -1.151  0.24981    
# dependency2:environment3        -1.201e-01  1.913e-02  4.876e+03  -6.279  3.7e-10 ***
# dependency2:group2               2.667e-02  2.212e-02  2.400e+02   1.205  0.22920    
# dependency2:group3               5.933e-02  2.191e-02  2.412e+02   2.708  0.00725 ** 
# environment2:group2             -1.757e-02  2.534e-02  3.433e+02  -0.693  0.48852    
# environment3:group2             -1.713e-02  2.850e-02  2.168e+02  -0.601  0.54840    
# environment2:group3             -3.563e-02  2.512e-02  3.460e+02  -1.418  0.15698    
# environment3:group3              7.733e-03  2.828e-02  2.197e+02   0.273  0.78478    
# dependency2:environment2:group2  7.756e-02  4.721e-02  4.826e+03   1.643  0.10044    
# dependency2:environment3:group2  7.115e-02  4.749e-02  4.864e+03   1.498  0.13410    
# dependency2:environment2:group3  9.442e-02  4.677e-02  4.860e+03   2.019  0.04357 *  
# dependency2:environment3:group3  1.295e-01  4.699e-02  4.875e+03   2.756  0.00588 ** 

# table ...

dt <- summary(model1)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_spr_question_logrt_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_question_logrt_model.csv')
dt <- read_rds('tables/src_spr_question_logrt_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                |Estimate|$SE$ |$t$    |$p$    |
# |:----------------------------------------|:-------|:----|:------|:------|
# |(Intercept)                              |8.47    |0.02 |458.38 |< .001 | ***
# |Dependency[RP]                           |−0.07   |0.01 |−7.93  |< .001 | ***
# |Environment[Long]                        |0.01    |0.01 |1.05   |0.295  |
# |Environment[Island]                      |0.11    |0.01 |7.94   |< .001 | ***
# |Group[KLE]                               |−0.10   |0.04 |−2.38  |0.018  | *
# |Group[MLE]                               |0.12    |0.04 |2.82   |0.005  | **
# |Dependency[RP]:Environment[Long]         |−0.02   |0.02 |−1.15  |0.250  |
# |Dependency[RP]:Environment[Island]       |−0.12   |0.02 |−6.28  |< .001 | ***
# |Dependency[RP]:Group[KLE]                |0.03    |0.02 |1.21   |0.229  |
# |Dependency[RP]:Group[MLE]                |0.06    |0.02 |2.71   |0.007  | **
# |Environment[Long]:Group[KLE]             |−0.02   |0.03 |−0.69  |0.489  |
# |Environment[Island]:Group[KLE]           |−0.02   |0.03 |−0.60  |0.548  |
# |Environment[Long]:Group[MLE]             |−0.04   |0.03 |−1.42  |0.157  |
# |Environment[Island]:Group[MLE]           |0.01    |0.03 |0.27   |0.785  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |0.08    |0.05 |1.64   |0.100  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.07    |0.05 |1.50   |0.134  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |0.09    |0.05 |2.02   |0.044  | *
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.13    |0.05 |2.76   |0.006  | **

# joint tests ...

joint_tests(model1)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  62.866  <.0001 ***
# environment                    2 Inf  35.178  <.0001 ***
# group                          2 Inf  14.421  <.0001 ***
# dependency:environment         2 Inf  22.244  <.0001 ***
# dependency:group               2 Inf   3.697  0.0248 *
# environment:group              4 Inf   0.998  0.4072
# dependency:environment:group   4 Inf   2.182  0.0683 .

# post-hoc tests ...

pairwise <- model1 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS    -0.0126 0.0253 Inf  -0.496  0.9998
# pronoun - gap long        ENS    -0.0918 0.0253 Inf  -3.637  0.0025 **
# pronoun - gap island      ENS    -0.1996 0.0255 Inf  -7.818  <.0001 ***
# pronoun - gap short       KLE    -0.0355 0.0244 Inf  -1.453  0.7569
# pronoun - gap long        KLE    -0.0372 0.0245 Inf  -1.517  0.7102
# pronoun - gap island      KLE    -0.1513 0.0247 Inf  -6.121  <.0001 ***
# pronoun - gap short       MLE    -0.0279 0.0238 Inf  -1.171  0.9157
# pronoun - gap long        MLE    -0.0127 0.0241 Inf  -0.527  0.9997
# pronoun - gap island      MLE    -0.0854 0.0243 Inf  -3.511  0.0040 **

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_question_logrt_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_question_logrt_posthoc.csv')
dt <- read_rds('tables/src_spr_question_logrt_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate|$SE$ |$z$   |$p$    |
# |:-----|:-----------|:--------|:-------|:----|:-----|:------|
# |ENS   |Short       |RP − Gap |−0.01   |0.03 |−0.50 |> .999 |
# |ENS   |Long        |RP − Gap |−0.09   |0.03 |−3.64 |0.002  | **
# |ENS   |Island      |RP − Gap |−0.20   |0.03 |−7.82 |< .001 | ***
# |KLE   |Short       |RP − Gap |−0.04   |0.02 |−1.45 |0.757  |
# |KLE   |Long        |RP − Gap |−0.04   |0.02 |−1.52 |0.710  |
# |KLE   |Island      |RP − Gap |−0.15   |0.02 |−6.12 |< .001 | ***
# |MLE   |Short       |RP − Gap |−0.03   |0.02 |−1.17 |0.916  |
# |MLE   |Long        |RP − Gap |−0.01   |0.02 |−0.53 |> .999 |
# |MLE   |Island      |RP − Gap |−0.09   |0.02 |−3.51 |0.004  | **

# ---------------------------------------------------------------------------- #
# model: question RT proficiency ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- spr_qrt %>%
  filter(group %in% c('KLE', 'MLE')) %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(dependency = fct_drop(dependency),
         environment = fct_drop(environment)) %>%
  mutate_at(c('study', 'group', 'participant', 'item'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

# check participants ...

check <- md %>%
  group_by(study, group) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant, proficiency, dependency, environment) %>%
  summarise(mean = mean(question_rt, na.rm = TRUE),
            mean2 = mean(log(question_rt), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

d1 <- aes(x = proficiency, y = mean, col = dependency, fill = dependency)
d2 <- aes(x = proficiency, y = mean2, col = dependency, fill = dependency)

s <- list(
  geom_smooth(aes(linetype = dependency), 
              method=lm, 
              alpha = .75,
              linewidth = .5),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        axis.title.x = element_markdown(),
        strip.text.y = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(group_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + d1 + s + 
  scale_y_continuous(name="RT (ms)")
ggsave('plots/orc_spr_question_rt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + d1 + s + 
  scale_y_continuous(name="RT (ms)")
ggsave('plots/src_spr_question_rt_proficiency.png', width=6.5, height=3, dpi=600)

p1 + d2 + s + 
  coord_cartesian(ylim = c(8, 9)) +
  scale_y_continuous(name="Log RT", breaks = c(8, 8.5, 9))
ggsave('plots/orc_spr_question_logrt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + d2 + s + 
  coord_cartesian(ylim = c(8, 9)) +
  scale_y_continuous(name="Log RT", breaks = c(8, 8.5, 9))
ggsave('plots/src_spr_question_logrt_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model2 <- lmer(log(question_rt) ~ dependency * environment * group * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 + dependency + environment | item),
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 33.83 sec elapsed

check_convergence(model2) # okay
check_singularity(model2) # warning (but random effects table looks okay)

# save ...

model2 %>% write_rds('models/orc_spr_question_logrt_proficiency_model.rds')
model2 <- read_rds('models/orc_spr_question_logrt_proficiency_model.rds')

# summary ...

summary(model2)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                  8.461e+00  1.992e-02  1.414e+02 424.718  < 2e-16 ***
# dependency2                                 -8.928e-02  1.318e-02  3.754e+01  -6.775 5.28e-08 ***
# environment2                                 3.769e-02  1.640e-02  3.222e+01   2.298  0.02817 *  
# environment3                                 3.135e-02  1.417e-02  3.680e+01   2.213  0.03316 *  
# group2                                       1.240e-01  3.686e-02  1.309e+02   3.364  0.00101 ** 
# proficiency                                 -5.869e-02  1.854e-02  1.315e+02  -3.166  0.00192 ** 
# dependency2:environment2                    -2.391e-02  2.395e-02  3.418e+03  -0.999  0.31807    
# dependency2:environment3                     7.003e-03  2.393e-02  3.386e+03   0.293  0.76981    
# dependency2:group2                           1.592e-03  2.107e-02  1.595e+02   0.076  0.93988    
# environment2:group2                          3.776e-02  2.657e-02  1.616e+02   1.421  0.15722    
# environment3:group2                          2.593e-02  2.606e-02  2.000e+02   0.995  0.32090    
# dependency2:proficiency                      1.708e-02  1.081e-02  1.712e+02   1.580  0.11588    
# environment2:proficiency                    -1.241e-02  1.348e-02  1.687e+02  -0.921  0.35861    
# environment3:proficiency                    -5.813e-03  1.324e-02  2.097e+02  -0.439  0.66120    
# group2:proficiency                           5.779e-02  3.707e-02  1.315e+02   1.559  0.12141    
# dependency2:environment2:group2             -5.204e-02  4.778e-02  3.408e+03  -1.089  0.27611    
# dependency2:environment3:group2             -4.507e-02  4.783e-02  3.412e+03  -0.942  0.34618    
# dependency2:environment2:proficiency         2.006e-02  2.453e-02  3.424e+03   0.818  0.41359    
# dependency2:environment3:proficiency         1.725e-02  2.462e-02  3.448e+03   0.701  0.48353    
# dependency2:group2:proficiency               1.457e-02  2.157e-02  1.703e+02   0.676  0.50025    
# environment2:group2:proficiency              3.245e-02  2.697e-02  1.689e+02   1.203  0.23050    
# environment3:group2:proficiency              5.192e-02  2.650e-02  2.099e+02   1.959  0.05144 .  
# dependency2:environment2:group2:proficiency -1.976e-02  4.891e-02  3.444e+03  -0.404  0.68629    
# dependency2:environment3:group2:proficiency  8.546e-02  4.908e-02  3.455e+03   1.741  0.08171 . 

# table ...

dt <- summary(model2)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_spr_question_logrt_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_question_logrt_proficiency_model.csv')
dt <- read_rds('tables/orc_spr_question_logrt_proficiency_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                            |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                                          |8.46                     |0.02 |424.72 |< .001 | ***
# |Dependency[RP]                                       |−0.09                    |0.01 |−6.77  |< .001 | ***
# |Environment[Long]                                    |0.04                     |0.02 |2.30   |0.028  | *
# |Environment[Island]                                  |0.03                     |0.01 |2.21   |0.033  | *
# |Group[KLE]                                           |0.12                     |0.04 |3.36   |0.001  | **
# |Proficiency                                          |−0.06                    |0.02 |−3.17  |0.002  | **
# |Dependency[RP]:Environment[Long]                     |−0.02                    |0.02 |−1.00  |0.318  |
# |Dependency[RP]:Environment[Island]                   |0.01                     |0.02 |0.29   |0.770  |
# |Dependency[RP]:Group[KLE]                            |0.00                     |0.02 |0.08   |0.940  |
# |Environment[Long]:Group[KLE]                         |0.04                     |0.03 |1.42   |0.157  |
# |Environment[Island]:Group[KLE]                       |0.03                     |0.03 |1.00   |0.321  |
# |Dependency[RP]:Proficiency                           |0.02                     |0.01 |1.58   |0.116  |
# |Environment[Long]:Proficiency                        |−0.01                    |0.01 |−0.92  |0.359  |
# |Environment[Island]:Proficiency                      |−0.01                    |0.01 |−0.44  |0.661  |
# |Group[KLE]:Proficiency                               |0.06                     |0.04 |1.56   |0.121  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]               |−0.05                    |0.05 |−1.09  |0.276  |
# |Dependency[RP]:Envir.[Island]:Group[KLE]             |−0.05                    |0.05 |−0.94  |0.346  |
# |Dependency[RP]:Envir.[Long]:Proficiency              |0.02                     |0.02 |0.82   |0.414  |
# |Dependency[RP]:Envir.[Island]:Proficiency            |0.02                     |0.02 |0.70   |0.484  |
# |Dependency[RP]:Group[KLE]:Proficiency                |0.01                     |0.02 |0.68   |0.500  |
# |Environment[Long]:Group[KLE]:Proficiency             |0.03                     |0.03 |1.20   |0.230  |
# |Environment[Island]:Group[KLE]:Proficiency           |0.05                     |0.03 |1.96   |0.051  | .
# |Dependency[RP]:Envir.[Long]:Group[KLE]:Proficiency   |−0.02                    |0.05 |−0.40  |0.686  |
# |Dependency[RP]:Envir.[Island]:Group[KLE]:Proficiency |0.09                     |0.05 |1.74   |0.082  | .

# joint tests ...

joint_tests(model2)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  45.718  <.0001 ***
# environment                                2 Inf   3.500  0.0302 *
# group                                      1 Inf  11.428  0.0007 ***
# proficiency                                1 Inf  10.027  0.0015 **
# dependency:environment                     2 Inf   0.908  0.4034
# dependency:group                           1 Inf   0.007  0.9340
# dependency:proficiency                     1 Inf   2.497  0.1140
# environment:group                          2 Inf   1.158  0.3142
# environment:proficiency                    2 Inf   0.428  0.6518
# group:proficiency                          1 Inf   2.430  0.1190
# dependency:environment:group               2 Inf   0.693  0.4999
# dependency:environment:proficiency         2 Inf   0.391  0.6766
# dependency:group:proficiency               1 Inf   0.456  0.4993
# environment:group:proficiency              2 Inf   2.038  0.1303
# dependency:environment:group:proficiency   2 Inf   2.604  0.0740 .

# post-hoc tests ...

pairwise <- model2 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group  estimate     SE  df z.ratio p.value
# pronoun - gap short       KLE    0.008311 0.0255 Inf   0.326  0.9997
# pronoun - gap long        KLE    0.038246 0.0252 Inf   1.521  0.5609
# pronoun - gap island      KLE   -0.017168 0.0255 Inf  -0.673  0.9844
# pronoun - gap short       MLE    0.000978 0.0249 Inf   0.039  1.0000
# pronoun - gap long        MLE    0.011157 0.0250 Inf   0.446  0.9983
# pronoun - gap island      MLE    0.060959 0.0250 Inf   2.440  0.0848 .

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_question_logrt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_question_logrt_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_spr_question_logrt_proficiency_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:------|
# |KLE   |Short       |RP − Gap |0.01                     |0.03 |0.33  |> .999 |
# |KLE   |Long        |RP − Gap |0.04                     |0.03 |1.52  |0.561  |
# |KLE   |Island      |RP − Gap |−0.02                    |0.03 |−0.67 |0.984  |
# |MLE   |Short       |RP − Gap |0.00                     |0.02 |0.04  |> .999 |
# |MLE   |Long        |RP − Gap |0.01                     |0.03 |0.45  |0.998  |
# |MLE   |Island      |RP − Gap |0.06                     |0.02 |2.44  |0.085  | *

# plot ...

emmip(model2, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model 5 (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model5 <- lmer(log(question_rt) ~ dependency * environment * group * proficiency + 
                 (1 + dependency * environment | participant) +
                 (1 | item),
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model5)
toc()
beep()

# 21.97 sec elapsed
# singular fit

check_convergence(model5) # warning (but lme4 says it's okay)
check_singularity(model5) # warning (but random effects table looks okay)

# save ...

model5 %>% write_rds('models/src_spr_question_logrt_proficiency_model.rds')
model5 <- read_rds('models/src_spr_question_logrt_proficiency_model.rds')

# summary ...

summary(model5)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                                   8.476575   0.022375 135.334606 378.844  < 2e-16 ***
# dependency2                                  -0.061489   0.011525 151.615666  -5.335 3.41e-07 ***
# environment2                                 -0.001475   0.013128 233.434264  -0.112  0.91066    
# environment3                                  0.108750   0.014143 154.346160   7.689 1.62e-12 ***
# group2                                        0.224716   0.042615 122.893649   5.273 5.83e-07 ***
# proficiency                                  -0.002779   0.021647 125.304325  -0.128  0.89804    
# dependency2:environment2                      0.021559   0.024557 676.878909   0.878  0.38030    
# dependency2:environment3                     -0.080907   0.027903 133.929639  -2.900  0.00437 ** 
# dependency2:group2                            0.035664   0.023083 152.231952   1.545  0.12441    
# environment2:group2                          -0.026179   0.026262 233.585226  -0.997  0.31988    
# environment3:group2                           0.017198   0.028292 154.405821   0.608  0.54417    
# dependency2:proficiency                       0.005510   0.012000 159.939648   0.459  0.64676    
# environment2:proficiency                     -0.009025   0.013654 248.370324  -0.661  0.50924    
# environment3:proficiency                     -0.014426   0.014769 161.831647  -0.977  0.33015    
# group2:proficiency                            0.070951   0.043292 125.296303   1.639  0.10375    
# dependency2:environment2:group2               0.003843   0.048971 668.593420   0.078  0.93748    
# dependency2:environment3:group2               0.081971   0.055858 134.470413   1.467  0.14458    
# dependency2:environment2:proficiency         -0.033279   0.025500 695.163688  -1.305  0.19230    
# dependency2:environment3:proficiency          0.034342   0.029072 140.614370   1.181  0.23948    
# dependency2:group2:proficiency               -0.022143   0.024000 159.916681  -0.923  0.35758    
# environment2:group2:proficiency              -0.010670   0.027430 251.209926  -0.389  0.69761    
# environment3:group2:proficiency               0.004828   0.029651 163.520307   0.163  0.87085    
# dependency2:environment2:group2:proficiency   0.093125   0.050997 695.173235   1.826  0.06827 .  
# dependency2:environment3:group2:proficiency   0.038167   0.058147 140.622096   0.656  0.51265   

# table ...

dt <- summary(model5)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_spr_question_logrt_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_question_logrt_proficiency_model.csv')
dt <- read_rds('tables/src_spr_question_logrt_proficiency_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                            |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                                          |8.48                     |0.02 |378.84 |< .001 | ***
# |Dependency[RP]                                       |−0.06                    |0.01 |−5.34  |< .001 | ***
# |Environment[Long]                                    |0.00                     |0.01 |−0.11  |0.911  |
# |Environment[Island]                                  |0.11                     |0.01 |7.69   |< .001 | ***
# |Group[KLE]                                           |0.22                     |0.04 |5.27   |< .001 | ***
# |Proficiency                                          |0.00                     |0.02 |−0.13  |0.898  |
# |Dependency[RP]:Environment[Long]                     |0.02                     |0.02 |0.88   |0.380  |
# |Dependency[RP]:Environment[Island]                   |−0.08                    |0.03 |−2.90  |0.004  | **
# |Dependency[RP]:Group[KLE]                            |0.04                     |0.02 |1.55   |0.124  |
# |Environment[Long]:Group[KLE]                         |−0.03                    |0.03 |−1.00  |0.320  |
# |Environment[Island]:Group[KLE]                       |0.02                     |0.03 |0.61   |0.544  |
# |Dependency[RP]:Proficiency                           |0.01                     |0.01 |0.46   |0.647  |
# |Environment[Long]:Proficiency                        |−0.01                    |0.01 |−0.66  |0.509  |
# |Environment[Island]:Proficiency                      |−0.01                    |0.01 |−0.98  |0.330  |
# |Group[KLE]:Proficiency                               |0.07                     |0.04 |1.64   |0.104  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]               |0.00                     |0.05 |0.08   |0.937  |
# |Dependency[RP]:Envir.[Island]:Group[KLE]             |0.08                     |0.06 |1.47   |0.145  |
# |Dependency[RP]:Envir.[Long]:Proficiency              |−0.03                    |0.03 |−1.31  |0.192  |
# |Dependency[RP]:Envir.[Island]:Proficiency            |0.03                     |0.03 |1.18   |0.239  |
# |Dependency[RP]:Group[KLE]:Proficiency                |−0.02                    |0.02 |−0.92  |0.358  |
# |Environment[Long]:Group[KLE]:Proficiency             |−0.01                    |0.03 |−0.39  |0.698  |
# |Environment[Island]:Group[KLE]:Proficiency           |0.00                     |0.03 |0.16   |0.871  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]:Proficiency   |0.09                     |0.05 |1.83   |0.068  | .
# |Dependency[RP]:Envir.[Island]:Group[KLE]:Proficiency |0.04                     |0.06 |0.66   |0.513  |

# joint tests ...

joint_tests(model5)

# model term                               df1 df2 F.ratio p.value
# dependency                                 1 Inf  28.506  <.0001 ***
# environment                                2 Inf  36.991  <.0001 ***
# group                                      1 Inf  27.648  <.0001 ***
# proficiency                                1 Inf   0.016  0.8978
# dependency:environment                     2 Inf   7.420  0.0006 ***
# dependency:group                           1 Inf   2.411  0.1205
# dependency:proficiency                     1 Inf   0.211  0.6461
# environment:group                          2 Inf   1.159  0.3138
# environment:proficiency                    2 Inf   0.512  0.5991
# group:proficiency                          1 Inf   2.686  0.1012
# dependency:environment:group               2 Inf   1.311  0.2694
# dependency:environment:proficiency         2 Inf   2.845  0.0581 .
# dependency:group:proficiency               1 Inf   0.851  0.3562
# environment:group:proficiency              2 Inf   0.143  0.8668
# dependency:environment:group:proficiency   2 Inf   1.687  0.1851

# post-hoc tests ...

pairwise <- model5 %>%
  emtrends(~ dependency * environment * group, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       KLE     0.0381 0.0259 Inf   1.473  0.5965
# pronoun - gap long        KLE    -0.0417 0.0269 Inf  -1.551  0.5375
# pronoun - gap island      KLE     0.0534 0.0327 Inf   1.634  0.4761
# pronoun - gap short       MLE    -0.0278 0.0244 Inf  -1.140  0.8278
# pronoun - gap long        MLE    -0.0145 0.0266 Inf  -0.546  0.9949
# pronoun - gap island      MLE     0.0256 0.0324 Inf   0.790  0.9654

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_question_logrt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_question_logrt_proficiency_posthoc.csv')
dt <- read_rds('tables/src_spr_question_logrt_proficiency_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$   |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:-----|
# |KLE   |Short       |RP − Gap |0.04                     |0.03 |1.47  |0.597 |
# |KLE   |Long        |RP − Gap |−0.04                    |0.03 |−1.55 |0.538 |
# |KLE   |Island      |RP − Gap |0.05                     |0.03 |1.63  |0.476 |
# |MLE   |Short       |RP − Gap |−0.03                    |0.02 |−1.14 |0.828 |
# |MLE   |Long        |RP − Gap |−0.01                    |0.03 |−0.55 |0.995 |
# |MLE   |Island      |RP − Gap |0.03                     |0.03 |0.79  |0.965 |

# plot ...

emmip(model5, dependency ~ proficiency | group + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: question RT (by accuracy) ####
# ---------------------------------------------------------------------------- #

# summarize to comprehension question data ...

temp1 <- spr_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  mutate(accuracy = as.logical(accuracy)) %>%
  group_by(study, group, group_label, participant, item, dependency, environment, accuracy, question_rt) %>%
  summarise() %>%
  ungroup()

# remove extreme values ...

temp2 <- temp1 %>%
  filter(question_rt <= 30000, question_rt >= 100)

# calculate percent excluded ...

check <- temp1 %>%
  mutate(keep = case_when(question_rt <= 30000 & question_rt >= 100 ~ TRUE,
                          TRUE ~ FALSE))

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# remove outliers more than 2 SDs from the mean ...

temp3 <- temp2 %>%
  group_by(study, group, dependency, environment, accuracy) %>%
  mutate(upper = mean(question_rt) + (2 * sd(question_rt)),
         lower = mean(question_rt) - (2 * sd(question_rt))) %>%
  ungroup() %>%
  filter(question_rt <= upper & question_rt >= lower) %>%
  select(-upper, -lower)

# calculate percent excluded ...

check <- temp2 %>%
  group_by(study, group, dependency, environment, accuracy) %>%
  mutate(upper = mean(question_rt) + (2 * sd(question_rt)),
         lower = mean(question_rt) - (2 * sd(question_rt))) %>%
  ungroup() %>%
  mutate(keep = case_when(question_rt <= upper & question_rt >= lower ~ TRUE,
                          TRUE ~ FALSE)) %>%
  select(study, group, participant, dependency, environment, accuracy, question_rt, upper, lower, keep)

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, group) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# plot (by accuracy) ...

plot <- temp3 %>%
  group_by(study, group, group_label, dependency, environment, accuracy) %>%
  summarise(mean = mean(question_rt, na.rm=T),
            lcl = my.lcl(question_rt),
            ucl = my.ucl(question_rt)) %>%
  ungroup()

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Environment", 
                   limits = c("short", "long", "island"), 
                   labels = c("Short", "Long", "Island")),
  scale_y_continuous(name="RT (ms)",
                     breaks = c(4000, 5000, 6000, 7000)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    strip.text = element_markdown(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin=margin(0, 0, 0, -5)),
  facet_grid2(accuracy~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(4000, 7000))
ggsave("plots/orc_spr_question_rt_accuracy.png", width=6.5, height=2.5, dpi=600)

p2 + s + coord_cartesian(ylim = c(4000, 7000))
ggsave("plots/src_spr_question_rt_accuracy.png", width=6.5, height=2.5, dpi=600)

# plot (accurate trials only) ...

plot <- plot %>%
  filter(accuracy == TRUE)

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Environment", 
                   limits = c("short", "long", "island"), 
                   labels = c("Short", "Long", "Island")),
  scale_y_continuous(name="RT (ms)",
                     breaks = c(4000, 5000, 6000, 7000)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    strip.text = element_markdown(),
    plot.margin = margin(2, 2, 2, 2),
    legend.margin=margin(0, 0, 0, -5)),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(4000, 6000))
ggsave("plots/orc_spr_question_rt_accurate.png", width=6.5, height=2.5, dpi=600)

p2 + s + coord_cartesian(ylim = c(4000, 7000))
ggsave("plots/src_spr_question_rt_accurate.png", width=6.5, height=2.5, dpi=600)

# ============================================================================ #
# ** ACCEPTABILITY JUDGMENT TASK ** ####
# ============================================================================ #

# create dataframe for ajt data ...

ajt <- df %>%
  filter(task %in% c('english_ajt', 'korean_ajt', 'mandarin_ajt')) %>%
  arrange(study, group, participant) %>%
  select_if(function(x){!all(is.na(x))})

# add pairing factor ...

ajt <- ajt %>%
  mutate(pairing = case_when(group == 'ENS' & task == 'english_ajt' ~ 'EE',
                         group == 'KLE' & task == 'english_ajt' ~ 'KE',
                         group == 'KLE' & task == 'korean_ajt' ~ 'KK',
                         group == 'MLE' & task == 'english_ajt' ~ 'ME',
                         group == 'MLE' & task == 'mandarin_ajt' ~ 'MM')) %>%
  mutate(pairing = fct_relevel(pairing, 'EE', 'KE', 'ME', 'KK', 'MM'))

# add language factor ...

ajt <- ajt %>%
  mutate(language = case_when(pairing %in% c('EE', 'KE', 'ME') ~ 'English',
                              pairing == 'KK' ~ 'Korean',
                              pairing == 'MM' ~ 'Mandarin'))

# check participants ...

check <- ajt %>%
  group_by(study, pairing) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

# check NA responses

check <- ajt %>%
  filter(is.na(response == TRUE)) %>%
  group_by(study, pairing) %>%
  summarise(n = n()) %>%
  ungroup()

# remove NA responses ...

ajt <- ajt %>%
  filter(is.na(response) == FALSE)

# add binary accept/reject column ...

ajt <- ajt %>%
  mutate(acceptance = case_when(response > 3.5 ~ TRUE,
                                response < 3.5 ~ FALSE))

# calculate z-scores for each task by participant on ratings for all trials ...

ajt <- ajt %>%
  mutate(response = as.numeric(response)) %>%
  group_by(study, group, task, participant) %>%
  mutate(zscore = (response - mean(response, na.rm=T)) / sd(response, na.rm = T)) %>%
  ungroup()

# label filler sets ...

ajt <- ajt %>%
  mutate(set = case_when(item %in% c('item31', 'item32', 'item33', 'item34', 'item35', 'item36', 'item37') ~ 'set1',
                         item %in% c('item38', 'item39', 'item40', 'item41', 'item42', 'item43', 'item44') ~ 'set2',
                         item %in% c('item45', 'item46', 'item47', 'item48', 'item49', 'item50', 'item51') ~ 'set3',
                         item %in% c('item52', 'item53', 'item54', 'item55', 'item56', 'item57', 'item58') ~ 'set4',
                         item %in% c('item59', 'item60', 'item61', 'item62', 'item63', 'item64', 'item65') ~ 'set5',
                         item %in% c('item66', 'item67', 'item68', 'item69', 'item70', 'item71', 'item72') ~ 'set6'))

# identify outliers ...

temp <- ajt %>%
  mutate(condition = case_when(pairing %in% c('EE', 'KE', 'ME') & dependency == 'pronoun' ~ 'ungrammatical',
                               pairing %in% c('EE', 'KE', 'ME') & environment == 'island' ~ 'ungrammatical',
                               pairing %in% c('EE', 'KE', 'ME') & dependency == 'gap' & environment %in% c('short', 'long') ~ 'grammatical',
                               pairing == 'KK' & dependency == 'pronoun' ~ 'ungrammatical',
                               pairing == 'KK' & dependency == 'gap' ~ 'grammatical',
                               pairing == 'MM' & dependency == 'pronoun' & environment == 'short' ~ 'ungrammatical',
                               pairing == 'MM' & dependency == 'pronoun' & environment %in% c('long', 'island') ~ 'ungrammatical',
                               pairing %in% c('KK', 'MM') & dependency == 'gap' ~ 'grammatical',
                               TRUE ~ condition)) %>%
  mutate_at(c('study', 'pairing', 'type', 'condition', 'participant'), factor) %>%
  group_by(study, pairing, participant, type, condition) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(study, pairing, type, condition) %>%
  mutate(mean = mean(response, na.rm = TRUE),
         sd = sd(response, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(outlier = case_when(condition == 'grammatical' & response < 3.5 & response < (mean - (3 * sd)) ~ TRUE,
                             condition == 'ungrammatical' & response > 3.5 & response > (mean + (3 * sd)) ~ TRUE)) %>%
  arrange(study, pairing, participant, type, condition) %>%
  group_by(participant) %>%
  mutate(exclude = case_when(sum(outlier[condition == 'grammatical']) == 2 ~ 'TRUE',
                             sum(outlier[condition == 'ungrammatical']) == 2 ~ 'TRUE')) %>%
  ungroup() %>%
  filter(exclude == TRUE, is.na(outlier) == FALSE)

# remove outliers ...

ajt_clean <- ajt %>%
  filter(!participant %in% temp$participant)

# check participants ...

check <- ajt_clean %>%
  group_by(study, pairing) %>%
  summarise(pairing_size = length(unique(participant))) %>%
  ungroup()

# ---------------------------------------------------------------------------- #
# plot: initial inspection by item & by participant ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

temp <- ajt_clean %>%
  mutate(item = str_remove(item, 'item')) %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '\n(n = ', pairing_size, ')')) %>%
  mutate(`Factor 1` = case_when(dependency == 'gap' ~ 'Gap',
                                dependency == 'pronoun' ~ 'RP',
                                condition == 'grammatical' ~ 'Grammatical',
                                condition == 'ungrammatical' ~ 'Ungrammatical'),
         `Factor 2` = case_when(environment == 'short' ~ 'Short',
                                environment == 'long' ~ 'Long',
                                environment == 'island' ~ 'Island',
                                # TRUE ~ 'Filler',
                                set %in% c('set1', 'set2') ~ 'Filler1',
                                set %in% c('set3', 'set4') ~ 'Filler2',
                                set == 'set5' ~ 'Filler3',
                                set == 'set6' ~ 'Filler4')) %>%
  mutate(factor1 = fct_relevel(`Factor 1`, 'Gap', 'RP', 'Grammatical', 'Ungrammatical'),
         factor2 = fct_relevel(`Factor 2`, 'Short', 'Long', 'Island', 'Filler1', 'Filler2', 'Filler3', 'Filler4'))

# summarise by item ...

plot1 <- temp %>%
  group_by(study, pairing, pairing_label, `Factor 1`, `Factor 2`, item) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(study, pairing) %>%
  # mutate(width = .1 * n()) %>%
  ungroup()

# summarise by participant ...

plot2 <- temp %>%
  group_by(study, pairing, pairing_label, `Factor 1`, `Factor 2`, participant) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(study, pairing) %>%
  # mutate(width = .1 * n()) %>%
  ungroup()

# plot ...

p1 <- ggplot(data = filter(plot1, study == 'ORC'))
p2 <- ggplot(data = filter(plot1, study == 'SRC'))
p3 <- ggplot(data = filter(plot2, study == 'ORC'))
p4 <- ggplot(data = filter(plot2, study == 'SRC'))
p5 <- ggplot(data = filter(temp, study == 'ORC'))
p6 <- ggplot(data = filter(temp, study == 'SRC'))

design <- '
ABD
#CE
'

s <- list(
  aes(x = `Factor 2`, y = response, fill = `Factor 1`, col = `Factor 1`),
  geom_hline(yintercept=3.5),
  #geom_violin(position = position_dodge(width = .5),
              #alpha = .5, 
              #col = NA
              #),
  geom_boxplot(outlier.size = .5, 
               width = .25,
               position = position_dodge(width = .5, preserve = 'single'), 
               # fill = NA
               ),
  #geom_text(size = 2.5, col = "black", position = position_jitter(seed=2)),
  theme_classic(),
  scale_x_discrete(name = 'Environment',
                   limits = c('Short', 'Long', 'Island', 'Filler1', 'Filler2', 'Filler3', 'Filler4'),
                   labels = c('Sh.', 'Lo.', 'Is.', 'F1', 'F2', 'F3', 'F4')),
  scale_y_continuous(name="Rating", limits=c(1, 6), breaks=c(1, 2, 3, 4, 5, 6)),
  scale_color_manual(name = 'Type',
                     values=c('#648fff', '#ffb000', '#40B0A6', 'coral'),
                     limits = c('Gap', 'RP', 'Grammatical', 'Ungrammatical')),
  scale_fill_manual(name='Type', 
                    values=c('#c1d2ff', '#ffdf99', '#b0e2de', '#ffccb9'),
                    limits = c('Gap', 'RP', 'Grammatical', 'Ungrammatical')),
  theme(text = element_text(size = 12), 
        legend.position = c(.15, .22),
        axis.text = element_text(size = 9),
        legend.margin = margin(-10, 0, 0, 0),
        axis.title.y = element_text(hjust = .85)),
  facet_manual(~pairing_label, design = design, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave("plots/orc_ajt_rating_item.png", width=6.5, height=4, dpi=600)

p2 + s
ggsave("plots/src_ajt_rating_item.png", width=6.5, height=4, dpi=600)

p3 + s
ggsave("plots/orc_ajt_rating_ppt.png", width=6.5, height=4, dpi=600)

p4 + s
ggsave("plots/src_ajt_rating_ppt.png", width=6.5, height=4, dpi=600)

p5 + s
ggsave("plots/orc_ajt_rating_all.png", width=6.5, height=4, dpi=600)

p6 + s
ggsave("plots/src_ajt_rating_all.png", width=6.5, height=4, dpi=600)

# prep dataframe ...

temp2 <- temp %>%
  mutate(factor2 = fct_recode(factor2, 
                              'Filler' = 'Filler1',
                              'Filler' = 'Filler2',
                              'Filler' = 'Filler3',
                              'Filler' = 'Filler4'))

plot3 <- temp2 %>%
  group_by(study, pairing, pairing_label, factor1, factor2, item) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup()

plot4 <- temp2 %>%
  group_by(study, pairing, pairing_label, factor1, factor2, participant) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup()

# plot ...

p3a <- ggplot(data = filter(plot3, study == 'ORC'))
p3b <- ggplot(data = filter(plot3, study == 'SRC'))
p4a <- ggplot(data = filter(plot4, study == 'ORC'))
p4b <- ggplot(data = filter(plot4, study == 'SRC'))

design <- '
ABD
#CE
'

s <- list(
  aes(x = factor2, y = response, fill = factor1, col = factor1),
  geom_hline(yintercept=3.5),
  geom_violin(position = position_dodge(width = .5),
              col = NA,
              width = 2),
  geom_boxplot(outlier.size = .5, 
               width = .1,
               fill = 'white',
               position = position_dodge(width = .5, preserve = 'single'), 
               # fill = NA
  ),
  theme_classic(),
  scale_x_discrete(name = 'Environment',
                   limits = c('Short', 'Long', 'Island', 'Filler')),
  scale_y_continuous(name="Rating", limits=c(1, 6), breaks=c(1, 3.5, 6)),
  scale_color_manual(name = 'Type',
                     values=c('#648fff', '#ffb000', '#40B0A6', 'coral'),
                     limits = c('Gap', 'RP', 'Grammatical', 'Ungrammatical')),
  scale_fill_manual(name='Type', 
                    values=c('#c1d2ff', '#ffdf99', '#b0e2de', '#ffccb9'),
                    limits = c('Gap', 'RP', 'Grammatical', 'Ungrammatical')),
  theme(text = element_text(size = 12), 
        plot.margin = margin(2, 2, 0, 2),
        legend.position = c(.15, .22),
        axis.text = element_text(size = 9),
        legend.margin = margin(-10, 0, 0, 0),
        axis.title.y = element_text(hjust = .85)),
  facet_manual(~pairing_label, design = design, axes = 'all', remove_labels = FALSE)
) 

p3a + s
ggsave("plots/orc_ajt_rating_item_v2.png", width=6.5, height=4, dpi=600)

p3b + s
ggsave("plots/src_ajt_rating_item_v2.png", width=6.5, height=4, dpi=600)

p4a + s
ggsave("plots/orc_ajt_rating_ppt_v2.png", width=6.5, height=4, dpi=600)

p4b + s
ggsave("plots/src_ajt_rating_ppt_v2.png", width=6.5, height=4, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: rating ####
# ---------------------------------------------------------------------------- #

# filter ... 

crit <- ajt_clean %>% 
  filter(!condition %in% c('grammatical', 'ungrammatical'))

# summarize ...

plot <- crit %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(', pairing, '; *n* = ', pairing_size, ')')) %>%
  mutate(group_label = paste0(group, ' (*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, group_label, task, dependency, environment) %>%
  summarise(mean = mean(response, na.rm = T),
            lcl = my.lcl(response),
            ucl = my.ucl(response)) %>%
  ungroup() %>%
  filter(dependency != 'NA') %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = case_when(environment == 'Island' & pairing %in% c('KK', 'MM') ~ "‘Island’",
                                 TRUE ~ environment)) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island", "‘Island’")) %>%
  arrange(study, pairing, dependency, environment)

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))

p2 <- ggplot(data=filter(plot, study == 'SRC'))

p1_bucld <- ggplot(data=filter(plot, study == 'ORC', pairing %in% c('EE', 'KE', 'ME')))

p2_bucld <- ggplot(data=filter(plot, study == 'SRC', pairing %in% c('EE', 'KE', 'ME')))

design <- '
 ABD
 #CE
'

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, 
      shape=dependency),
  geom_hline(yintercept = 3.5),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                width = 0.6, 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name='Environment'),
  scale_y_continuous(name='Rating', 
                     limits=c(1, 6), 
                     breaks = c(1, 3.5, 6)),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency', 
                     values=c('41', 'solid'), 
                     labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
  #theme(
        #text = element_text(size = 12),
        #plot.margin = margin(2, 2, 2, 2),
        #legend.position = c(.15, .22),
        #legend.margin = margin(-10, 0, 0, 0),
        #axis.title.y = element_text(hjust = .88),
        #strip.text = element_markdown()),
  #facet_manual(~pairing_label, design = 'ABDCE', axes = 'all', 
               #remove_labels = FALSE, scales = 'free')
) 

p1 + s
ggsave("plots/orc_ajt_rating.png", width=6.5, height=2.25, dpi=600)

p2 + s
ggsave("plots/src_ajt_rating.png", width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot %>%
  select(study, pairing, dependency, environment, mean, lcl, ucl) %>%
  arrange(study, pairing, dependency, environment) %>%
  my_table() 

dt %>% write_rds('tables/ajt_rating.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/ajt_rating.csv')
dt <- read_rds('tables/ajt_rating.rds')

# ---------------------------------------------------------------------------- #
# model: rating ####
# ---------------------------------------------------------------------------- #

# prep dataframe

md <- ajt_clean %>%
  filter(!condition %in% c('grammatical', 'ungrammatical')) %>%
  select(study, pairing, participant, item, dependency, environment, response) %>%
  mutate_at(c('study', 'pairing', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(response = ordered(response))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol = 2)
contrasts(md$pairing) <- contr.treatment(5) - matrix(rep(1/5, 20), ncol = 4)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- clmm(response ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + environment + pairing | item), 
               data = filter(md, study == 'ORC'), control = clmm.control(grtol = 1e6))
summary(model7)
toc()
beep()

performance::check_singularity(model7) # FALSE
ordinal::convergence(model7)

broom::tidy(model7, conf.int = TRUE)

# 23814.66
# no warnings

model7 %>% write_rds('models/orc_ajt_rating_model')
model7 <- read_rds('models/orc_ajt_rating_model')

# model summary ...

summary(model7)

# Parameter                         Estimate Std. Error z value Pr(>|z|)    
# dependency2                       -2.28300    0.15254 -14.967  < 2e-16 ***
# environment2                      -0.33376    0.06970  -4.788 1.68e-06 ***
# environment3                      -0.53842    0.09828  -5.479 4.29e-08 ***
# pairing2                           0.37326    0.20553   1.816 0.069351 .  
# pairing3                           0.69318    0.20156   3.439 0.000584 ***
# pairing4                           0.25829    0.20564   1.256 0.209107    
# pairing5                          -0.39099    0.20110  -1.944 0.051861 .  
# dependency2:environment2           1.39134    0.10021  13.885  < 2e-16 ***
# dependency2:environment3           2.63157    0.10123  25.995  < 2e-16 ***
# dependency2:pairing2               1.46846    0.35451   4.142 3.44e-05 ***
# dependency2:pairing3               2.15675    0.34443   6.262 3.80e-10 ***
# dependency2:pairing4               1.61127    0.35438   4.547 5.45e-06 ***
# dependency2:pairing5               2.96488    0.34413   8.615  < 2e-16 ***
# environment2:pairing2              0.53364    0.18337   2.910 0.003613 ** 
# environment3:pairing2              1.91151    0.24386   7.838 4.56e-15 ***
# environment2:pairing3              0.32530    0.17807   1.827 0.067728 .  
# environment3:pairing3              1.15608    0.23584   4.902 9.48e-07 ***
# environment2:pairing4             -0.19341    0.18609  -1.039 0.298663    
# environment3:pairing4              1.10166    0.24514   4.494 6.99e-06 ***
# environment2:pairing5              0.77692    0.17682   4.394 1.11e-05 ***
# environment3:pairing5              0.99378    0.23482   4.232 2.32e-05 ***
# dependency2:environment2:pairing2 -0.44153    0.31903  -1.384 0.166358    
# dependency2:environment3:pairing2 -2.47640    0.31868  -7.771 7.80e-15 ***
# dependency2:environment2:pairing3 -0.67364    0.30858  -2.183 0.029030 *  
# dependency2:environment3:pairing3 -3.46304    0.30822 -11.235  < 2e-16 ***
# dependency2:environment2:pairing4  2.24384    0.32484   6.907 4.93e-12 ***
# dependency2:environment3:pairing4 -0.02598    0.32447  -0.080 0.936188    
# dependency2:environment2:pairing5 -0.43992    0.30381  -1.448 0.147619    
# dependency2:environment3:pairing5 -3.02971    0.30266 -10.010  < 2e-16 ***

# table ...

dt <- summary(model7)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_ajt_rating_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rating_model.csv')
dt <- read_rds('tables/orc_ajt_rating_model.rds')

# joint tests ...

joint_tests(model7)

# model term                     df1 df2 F.ratio p.value
# dependency                       1 Inf 224.012  <.0001 ***
# environment                      2 Inf  17.646  <.0001 ***
# pairing                          4 Inf  46.961  <.0001 ***
# dependency:environment           2 Inf 338.281  <.0001 ***
# dependency:pairing               4 Inf  26.436  <.0001 ***
# environment:pairing              8 Inf  13.973  <.0001 ***
# dependency:environment:pairing   8 Inf  37.906  <.0001 ***

# post-hoc tests ...

pairwise <- model7 %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate    SE  df z.ratio p.value
# pronoun - gap short       EE        -5.818 0.289 Inf -20.153  <.0001 ***
# pronoun - gap long        EE        -4.564 0.270 Inf -16.886  <.0001 ***
# pronoun - gap island      EE        -1.387 0.260 Inf  -5.341  <.0001 ***
# pronoun - gap short       KE        -3.377 0.298 Inf -11.323  <.0001 ***
# pronoun - gap long        KE        -2.565 0.294 Inf  -8.731  <.0001 ***
# pronoun - gap island      KE        -1.423 0.295 Inf  -4.827  <.0001 ***
# pronoun - gap short       ME        -2.282 0.281 Inf  -8.134  <.0001 ***
# pronoun - gap long        ME        -1.702 0.278 Inf  -6.131  <.0001 ***
# pronoun - gap island      ME        -1.315 0.278 Inf  -4.729  <.0001 ***
# pronoun - gap short       KK        -4.946 0.311 Inf -15.928  <.0001 ***
# pronoun - gap long        KK        -1.449 0.287 Inf  -5.042  <.0001 ***
# pronoun - gap island      KK        -0.541 0.289 Inf  -1.876  0.4284
# pronoun - gap short       MM        -1.697 0.277 Inf  -6.118  <.0001 ***
# pronoun - gap long        MM        -0.883 0.275 Inf  -3.216  0.0149 *
# pronoun - gap island      MM        -0.296 0.275 Inf  -1.076  0.9414

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_rating_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rating_posthoc.csv')
dt <- read_rds('tables/orc_ajt_rating_posthoc.rds')

# plot ...

plot(emmeans(model7, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- clmm(response ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + environment + pairing | item), 
               data = filter(md, study == 'SRC'), control = clmm.control(grtol = 1e6))
summary(model7)
toc()
beep()

performance::check_singularity(model7) # FALSE

# 14100.47 sec elapsed (3.92 hours)
# no warnings

model7 %>% write_rds('models/src_ajt_rating_model7.rds')
model7 <- read_rds('models/src_ajt_rating_model7.rds')

# model summary ...

summary(model7)

# Parameter                         Estimate Std. Error z value Pr(>|z|)    
# dependency2                       -1.59300    0.10324 -15.431  < 2e-16 ***
# environment2                      -0.63612    0.10323  -6.162 7.18e-10 ***
# environment3                      -1.48440    0.12069 -12.299  < 2e-16 ***
# pairing2                           0.59757    0.21150   2.825 0.004723 ** 
# pairing3                           1.24213    0.21205   5.858 4.69e-09 ***
# pairing4                           1.09862    0.21404   5.133 2.86e-07 ***
# pairing5                           1.33189    0.21290   6.256 3.95e-10 ***
# dependency2:environment2           3.23076    0.11651  27.730  < 2e-16 ***
# dependency2:environment3           5.34182    0.12266  43.548  < 2e-16 ***
# dependency2:pairing2               2.05731    0.26498   7.764 8.22e-15 ***
# dependency2:pairing3               3.90111    0.26117  14.937  < 2e-16 ***
# dependency2:pairing4               1.04351    0.26632   3.918 8.92e-05 ***
# dependency2:pairing5               2.48781    0.25863   9.619  < 2e-16 ***
# environment2:pairing2              0.15022    0.27616   0.544 0.586456    
# environment3:pairing2              1.77340    0.30586   5.798 6.71e-09 ***
# environment2:pairing3             -0.10667    0.27140  -0.393 0.694306    
# environment3:pairing3              1.12786    0.30018   3.757 0.000172 ***
# environment2:pairing4              0.56947    0.27930   2.039 0.041454 *  
# environment3:pairing4              2.47755    0.30899   8.018 1.07e-15 ***
# environment2:pairing5              0.41335    0.26920   1.536 0.124661    
# environment3:pairing5              1.84157    0.29772   6.186 6.18e-10 ***
# dependency2:environment2:pairing2  0.63868    0.39244   1.627 0.103643    
# dependency2:environment3:pairing2 -2.33295    0.39672  -5.881 4.09e-09 ***
# dependency2:environment2:pairing3  0.81689    0.38564   2.118 0.034154 *  
# dependency2:environment3:pairing3 -3.77120    0.39160  -9.630  < 2e-16 ***
# dependency2:environment2:pairing4  2.26453    0.40359   5.611 2.01e-08 ***
# dependency2:environment3:pairing4 -3.25153    0.40773  -7.975 1.53e-15 ***
# dependency2:environment2:pairing5 -0.04792    0.38132  -0.126 0.899987    
# dependency2:environment3:pairing5 -4.58982    0.38600 -11.891  < 2e-16 ***

# table ...

dt <- summary(model7)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_ajt_rating_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rating_model.csv')
dt <- read_rds('tables/src_ajt_rating_model.rds')

# joint tests ...

joint_tests(model7)

# model term                     df1 df2 F.ratio p.value
# dependency                       1 Inf 238.105  <.0001 ***
# environment                      2 Inf  81.965  <.0001 ***
# pairing                          4 Inf  23.385  <.0001 ***
# dependency:environment           2 Inf 948.560  <.0001 ***
# dependency:pairing               4 Inf  87.001  <.0001 ***
# environment:pairing              8 Inf  15.031  <.0001 ***
# dependency:environment:pairing   8 Inf  47.107  <.0001 ***

# post-hoc tests ...

pairwise <- model7 %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate    SE  df z.ratio p.value
# pronoun - gap short       EE        -7.033 0.312 Inf -22.535  <.0001 ***
# pronoun - gap long        EE        -4.537 0.247 Inf -18.371  <.0001 ***
# pronoun - gap island      EE         1.098 0.246 Inf   4.456  0.0001 ***
# pronoun - gap short       KE        -4.411 0.234 Inf -18.891  <.0001 ***
# pronoun - gap long        KE        -1.276 0.215 Inf  -5.926  <.0001 ***
# pronoun - gap island      KE         1.387 0.215 Inf   6.439  <.0001 ***
# pronoun - gap short       ME        -2.147 0.215 Inf  -9.984  <.0001 ***
# pronoun - gap long        ME         1.166 0.205 Inf   5.674  <.0001 ***
# pronoun - gap island      ME         2.212 0.209 Inf  10.596  <.0001 ***
# pronoun - gap short       KK        -5.661 0.259 Inf -21.884  <.0001 ***
# pronoun - gap long        KK        -0.900 0.208 Inf  -4.324  0.0002 ***
# pronoun - gap island      KK        -0.781 0.209 Inf  -3.731  0.0027 **
# pronoun - gap short       MM        -3.000 0.214 Inf -13.994  <.0001 ***
# pronoun - gap long        MM        -0.551 0.201 Inf  -2.736  0.0798 .
# pronoun - gap island      MM         0.541 0.199 Inf   2.723  0.0826 .

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_rating_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rating_posthoc.csv')
dt <- read_rds('tables/src_ajt_rating_posthoc.rds')

# plot ...

plot(emmeans(model7, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: rating proficiency ####
# ---------------------------------------------------------------------------- #

# prep dataframe

md <- ajt_clean %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(!condition %in% c('grammatical', 'ungrammatical')) %>%
  select(study, pairing, participant, item, dependency, environment, response, proficiency) %>%
  mutate_at(c('study', 'pairing', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(response = ordered(response))

# remove participants with no proficiency scores ...

md <- md %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(proficiency = as.numeric(proficiency))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol = 2)
contrasts(md$pairing) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# plot ...

plot <- md %>%
  mutate(group = case_when(pairing == 'KE' ~ 'KLE', 
                           pairing == 'ME' ~ 'MLE')) %>%
  group_by(study, pairing) %>%
  mutate(response = as.numeric(as.character(response))) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, '<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, participant, proficiency, dependency, environment) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = response, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='Rating', breaks = c(1, 3.5, 6)),
  coord_cartesian(ylim = c(0.9, 6.1)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                    values=c('41', 'solid'), 
                    labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ajt_rating_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_ajt_rating_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #


tic()
model3 <- clmm(response ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency * environment | participant) + 
                 (0 + dependency * environment | item),
               data = filter(md, study == 'ORC'), control = clmm.control(grtol = 1e6))
summary(model3)
toc()
beep()

performance::check_singularity(model3) # FALSE

# 17497.84 sec elapsed (4.86 hours)
# no warnings

# save ...

model3 %>% write_rds('models/orc_ajt_rating_proficiency_model3.rds')
model3 <- read_rds('models/orc_ajt_rating_proficiency_model3.rds')

# model summary ...

summary(model3)

# Parameter                                     Estimate Std. Error z value Pr(>|z|)    
# proficiency                                    0.16929    0.13706   1.235 0.216793    
# dependency2                                   -2.84966    0.35970  -7.922 2.33e-15 ***
# environment2                                  -0.37623    0.17275  -2.178 0.029418 *  
# environment3                                  -0.20367    0.17793  -1.145 0.252348    
# pairing2                                       0.46339    0.27019   1.715 0.086332 .  
# proficiency:dependency2                       -1.41239    0.35594  -3.968 7.24e-05 ***
# proficiency:environment2                      -0.17228    0.15516  -1.110 0.266861    
# proficiency:environment3                      -0.47916    0.15709  -3.050 0.002286 ** 
# dependency2:environment2                       1.21267    0.27332   4.437 9.13e-06 ***
# dependency2:environment3                       2.50140    0.35564   7.033 2.01e-12 ***
# proficiency:pairing2                          -0.13970    0.27248  -0.513 0.608173    
# dependency2:pairing2                           1.10555    0.70173   1.575 0.115150    
# environment2:pairing2                         -0.18877    0.29730  -0.635 0.525464    
# environment3:pairing2                         -0.70281    0.30196  -2.327 0.019940 *  
# proficiency:dependency2:environment2           0.42572    0.22623   1.882 0.059863 .  
# proficiency:dependency2:environment3           1.35991    0.30981   4.390 1.14e-05 ***
# proficiency:dependency2:pairing2               1.24992    0.70840   1.764 0.077660 .  
# proficiency:environment2:pairing2             -0.13514    0.30177  -0.448 0.654295    
# proficiency:environment3:pairing2              0.08781    0.30686   0.286 0.774754    
# dependency2:environment2:pairing2             -0.61087    0.41869  -1.459 0.144563    
# dependency2:environment3:pairing2             -1.97231    0.59671  -3.305 0.000949 ***
# proficiency:dependency2:environment2:pairing2 -0.41745    0.43296  -0.964 0.334955    
# proficiency:dependency2:environment3:pairing2 -1.95664    0.60504  -3.234 0.001221 ** 

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_rating_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rating_proficiency_model.csv')
dt <- read_rds('tables/orc_ajt_rating_proficiency_model.rds')

# joint tests ...

joint_tests(model3)

# model term                                 df1 df2 F.ratio p.value
# proficiency                                  1 Inf   1.525  0.2168
# dependency                                   1 Inf  62.924  <.0001 ***
# environment                                  2 Inf   2.472  0.0844 .
# pairing                                      1 Inf   2.936  0.0866 .
# proficiency:dependency                       1 Inf  15.746  0.0001 ***
# proficiency:environment                      2 Inf   5.287  0.0051 **
# proficiency:pairing                          1 Inf   0.263  0.6082
# dependency:environment                       2 Inf  24.874  <.0001 ***
# dependency:pairing                           1 Inf   2.498  0.1140
# environment:pairing                          2 Inf   3.307  0.0366 *
# proficiency:dependency:environment           2 Inf  10.595  <.0001 ***
# proficiency:dependency:pairing               1 Inf   3.113  0.0777 .
# proficiency:environment:pairing              2 Inf   0.380  0.6838
# dependency:environment:pairing               2 Inf   5.959  0.0026 **
# proficiency:dependency:environment:pairing   2 Inf   6.442  0.0016 **

pairwise <- model3 %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment pairing estimate    SE  df z.ratio p.value
# pronoun - gap short       KE        -3.028 0.645 Inf  -4.697  <.0001 ***
# pronoun - gap long        KE        -2.394 0.562 Inf  -4.259  0.0001 ***
# pronoun - gap island      KE        -0.690 0.548 Inf  -1.258  0.5884
# pronoun - gap short       ME        -0.987 0.534 Inf  -1.848  0.2339
# pronoun - gap long        ME        -0.770 0.478 Inf  -1.612  0.3563
# pronoun - gap island      ME        -0.605 0.478 Inf  -1.268  0.5816

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_rating_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rating_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_ajt_rating_proficiency_posthoc.rds')

# plot ...

emmip(model3, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- clmm(response ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency * environment | participant) + 
                 (1 + dependency * environment | item), 
               data = filter(md, study == 'SRC'), control = clmm.control(grtol = 1e6))
summary(model7)
toc()
beep()

performance::check_singularity(model7) # FALSE

# no warnings

# save ...

model7 %>% write_rds('models/src_ajt_rating_proficiency_model7.rds')
model7 <- read_rds('models/src_ajt_rating_proficiency_model7.rds')

# summary ...

summary(model7)

# Parameter                                     Estimate Std. Error z value Pr(>|z|)    
# dependency2                                   -0.49457    0.23177  -2.134  0.03285 *  
# environment2                                  -1.10378    0.21852  -5.051 4.39e-07 ***
# environment3                                  -1.88401    0.22675  -8.309  < 2e-16 ***
# pairing2                                       0.87325    0.27962   3.123  0.00179 ** 
# proficiency                                    0.06103    0.14029   0.435  0.66354    
# dependency2:environment2                       4.69360    0.37503  12.515  < 2e-16 ***
# dependency2:environment3                       7.26396    0.51118  14.210  < 2e-16 ***
# dependency2:pairing2                           2.38312    0.44197   5.392 6.97e-08 ***
# environment2:pairing2                         -0.17831    0.40733  -0.438  0.66157    
# environment3:pairing2                         -0.70102    0.42520  -1.649  0.09922 .  
# dependency2:proficiency                       -0.44248    0.22120  -2.000  0.04547 *  
# environment2:proficiency                       0.33980    0.20479   1.659  0.09707 .  
# environment3:proficiency                       0.16957    0.21329   0.795  0.42660    
# pairing2:proficiency                          -0.03987    0.28016  -0.142  0.88683    
# dependency2:environment2:pairing2              0.25493    0.67961   0.375  0.70757    
# dependency2:environment3:pairing2             -1.25992    0.95063  -1.325  0.18505    
# dependency2:environment2:proficiency           0.75675    0.34071   2.221  0.02635 *  
# dependency2:environment3:proficiency           2.20211    0.47705   4.616 3.91e-06 ***
# dependency2:pairing2:proficiency               0.82297    0.44130   1.865  0.06220 .  
# environment2:pairing2:proficiency             -0.56274    0.40710  -1.382  0.16687    
# environment3:pairing2:proficiency             -0.40043    0.42409  -0.944  0.34507    
# dependency2:environment2:pairing2:proficiency  0.36758    0.67566   0.544  0.58642    
# dependency2:environment3:pairing2:proficiency -0.71479    0.94860  -0.754  0.45114  

# table ...

dt <- summary(model7)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_rating_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rating_proficiency_model.csv')
dt <- read_rds('tables/src_ajt_rating_proficiency_model.rds')

# joint tests ...

joint_tests(model7)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 Inf   4.559  0.0328 *
# environment                                  2 Inf  35.789  <.0001 ***
# pairing                                      1 Inf   9.753  0.0018 **
# proficiency                                  1 Inf   0.189  0.6635
# dependency:environment                       2 Inf 109.960  <.0001 ***
# dependency:pairing                           1 Inf  29.086  <.0001 ***
# dependency:proficiency                       1 Inf   4.001  0.0455 *
# environment:pairing                          2 Inf   1.928  0.1454
# environment:proficiency                      2 Inf   1.548  0.2126
# pairing:proficiency                          1 Inf   0.020  0.8868
# dependency:environment:pairing               2 Inf   1.975  0.1387
# dependency:environment:proficiency           2 Inf  10.937  <.0001 ***
# dependency:pairing:proficiency               1 Inf   3.478  0.0622 .
# environment:pairing:proficiency              2 Inf   0.958  0.3836
# dependency:environment:pairing:proficiency   2 Inf   1.077  0.3405

# post-hoc tests ...

pairwise <- model7 %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate    SE  df z.ratio p.value
# pronoun - gap short       KE       -1.8981 0.502 Inf  -3.778  0.0009 **
# pronoun - gap long        KE       -1.3252 0.488 Inf  -2.713  0.0382 *
# pronoun - gap island      KE        0.6614 0.433 Inf   1.527  0.5275
# pronoun - gap short       ME       -0.9594 0.445 Inf  -2.158  0.1634
# pronoun - gap long        ME       -0.0189 0.446 Inf  -0.042  1.0000
# pronoun - gap island      ME        0.8853 0.399 Inf   2.218  0.1421

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

table %>% write_rds('tables/src_ajt_rating_proficiency_posthoc.rds')
csv <- table %>% csv_table() %>% write_csv('tables/src_ajt_rating_proficiency_posthoc.csv')
table <- read_rds('tables/src_ajt_rating_proficiency_posthoc.rds')

# plot ...

emmip(model7, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: rating proficiency (ENSs included) ####
# ---------------------------------------------------------------------------- #

# prep dataframe

md <- ajt_clean %>%
  filter(pairing %in% c('EE', 'KE', 'ME')) %>%
  filter(!condition %in% c('grammatical', 'ungrammatical')) %>%
  select(study, pairing, participant, item, dependency, environment, response, proficiency) %>%
  mutate_at(c('study', 'pairing', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(response = ordered(response))

# remove participants with no proficiency scores ...

md <- md %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(proficiency = as.numeric(proficiency))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol = 2)
contrasts(md$pairing) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol = 2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# plot ...

plot <- md %>%
  mutate(group = case_when(pairing == 'EE' ~ 'ENS',
                           pairing == 'KE' ~ 'KLE', 
                           pairing == 'ME' ~ 'MLE')) %>%
  group_by(study, pairing) %>%
  mutate(response = as.numeric(as.character(response))) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, '<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, participant, proficiency, dependency, environment) %>%
  summarise(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = response, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='Rating', breaks = c(1, 3.5, 6)),
  coord_cartesian(ylim = c(0.9, 6.1)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ajt_rating_proficiency_ens.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_ajt_rating_proficiency_ens.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: rating filler ####
# ---------------------------------------------------------------------------- #

# filter to filler trials ...

fill <- ajt_clean %>% filter(condition %in% c('grammatical', 'ungrammatical'))

# create 'superset' factor ...

fill <- fill %>%
  mutate(superset = case_when(set %in% c('set1', 'set2') ~ 'superset1',
                              set %in% c('set3', 'set4') ~ 'superset2',
                              set == 'set5' ~ 'superset3',
                              set == 'set6' ~ 'superset4'))

# summarize for plotting ...

plot <- fill %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(', pairing, '; *n* = ', pairing_size, ')')) %>%
  mutate(condition = fct_drop(condition),
         superset = fct_drop(superset)) %>%
  group_by(study, pairing, pairing_label, set, superset, condition) %>%
  summarise(mean = mean(response, na.rm=T),
            lcl = my.lcl(response),
            ucl = my.ucl(response)) %>%
  mutate(width = 0.3 * n()) %>%
  ungroup() %>%
  filter(condition != 'NA')

# generate plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

design <- '
 ABD
 #CE
'

s <- list(
  aes(x=superset, y=mean, group=condition, col=condition, shape=condition, width = width * .75),
  geom_hline(yintercept = 3.5),
  geom_line(aes(linetype = condition), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                linetype = 1,
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Filler Type", 
                   limits = c("superset1", "superset2", "superset3", "superset4"), 
                   labels = c("F1", "F2", "F3", "F4")),
  scale_y_continuous(name='Rating', 
                     limits=c(1, 6), 
                     breaks = c(1, 3.5, 6)),
  scale_colour_manual(name='Grammaticality:', 
                      values=c('#40B0A6', 'coral'), 
                      labels=c('Grammatical', 'Ungrammatical')),
  scale_linetype_manual(name='Grammaticality:', 
                     values=c('solid', '41'), 
                     labels=c('Grammatical', 'Ungrammatical')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  # theme(text = element_text(size = 12),
  #       axis.text.x = element_markdown(),
  #       strip.text = element_markdown(),
  #       legend.position = c(.15, .22),
  #       legend.margin = margin(-10, 0, 0, 0),
  #       axis.title.y = element_text(hjust = .88),
  #       plot.title = element_text(size = 12, hjust = .5)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
) 

p1 + s
ggsave("plots/orc_ajt_filler_rating.png", width=6.5, height=2.25, dpi=600)

p2 + s
ggsave("plots/src_ajt_filler_rating.png", width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot %>%
  select(study, pairing, Type = set, Grammaticality = condition, mean, lcl, ucl) %>%
  my_table()

dt <- dt %>%
  mutate(Type = case_when(Type == "F1" ~ "F1G",
                          Type == "F2" ~ "F1U",
                          Type == "F3" ~ "F2G",
                          Type == "F4" ~ "F2U",
                          Type == "F5" ~ "F3G",
                          Type == "F6" ~ "F4G"))

dt %>% write_rds('tables/ajt_filler_rating.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/ajt_filler_rating.csv')
dt <- read_rds('tables/ajt_filler_rating.rds')

# ---------------------------------------------------------------------------- #
# plot: RP acceptors, rejectors, and inconsistent raters by environment ####
# ---------------------------------------------------------------------------- #

# categorize participants based on consistency of RP acceptance ... 

ajt_ind <- ajt_clean %>%
  filter(!condition %in% c('grammatical', 'ungrammatical')) %>%
  mutate(response = as.numeric(response))

ajt_ind <- ajt_ind %>%
  filter(dependency == 'pronoun') %>%
  group_by(study, group, language, pairing, participant, environment, proficiency) %>%
  summarise(acceptance = mean(acceptance, na.rm = TRUE),
            rating = mean(response, na.rm = TRUE)) %>%
  ungroup()

ajt_ind <- ajt_ind %>%
  mutate(type = case_when(acceptance >= .8 ~ 'accept',
                          acceptance <= .2 ~ 'reject',
                          TRUE ~ 'undefined')) %>%
  mutate(pairing = fct_relevel(pairing, 'EE', 'KE', 'ME', 'KK', 'MM')) %>%
  select(study, group, language, pairing, participant, environment, type, acceptance, rating, proficiency)

# plot ...

plot <- ajt_ind %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(', pairing, '; *n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, environment) %>%
  summarise(accept = length(which(type == 'accept')),
            reject = length(which(type == 'reject')),
            other = length(which(type == 'undefined')),
            total = n()) %>%
  select(study, pairing, pairing_label, environment, accept, reject, other, total) %>%
  ungroup() %>%
  pivot_longer(cols = c(accept:other), names_to = 'type', values_to = 'count') %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = case_when(environment == 'Island' & pairing %in% c('KK', 'MM') ~ "‘Island’",
                                 TRUE ~ environment)) %>%
  mutate(environment = fct_relevel(environment, 'Short', 'Long', 'Island', '‘Island’'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = environment, y = (count/total)*100, 
      group = type, col = type, 
      linetype = type, label = round(count/total * 100)),
  geom_hline(yintercept = 50),
  geom_line(position = position_dodge(width = .2)),
  geom_point(size = 1.5,
             position = position_dodge(width = .2), 
             shape = 21, 
             fill = 'white',
             show.legend = FALSE),
  # geom_text(size = 3, col = "black", hjust = 1.25, vjust = .5,
  #           position = position_dodge(width = .2)),
  theme_classic(),
  coord_cartesian(ylim = c(0, 100)),
  scale_x_discrete(name = 'Environment'),
  scale_y_continuous(name='Participants (%)', breaks = c(0, 50, 100)),
  scale_colour_manual(name='Category:', 
                      values=c('#40B0A6', 'coral', '#47b0e5'), 
                      limits = c('accept', 'reject', 'other'),
                      labels = c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  scale_linetype_manual(name='Category:', 
                        values=c('solid', '41', '5222'), 
                        limits = c('accept', 'reject', 'other'),
                        labels = c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(2, 2, 2, 2)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
)

p1 + s
ggsave('plots/orc_ajt_acceptor.png', width=6.5, height=2.25, dpi=600)

p2 + s
ggsave('plots/src_ajt_acceptor.png', width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot %>%
  mutate(Percent = count / total * 100) %>%
  select(study, pairing, environment, type, count, Percent) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'),
         type = fct_relevel(type, 'accept', 'reject', 'other')) %>%
  arrange(study, pairing, environment, type) %>%
  mutate_at(c('study', 'pairing', 'environment', 'type'), as.character) %>%
  mutate_all(str_replace_all, c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
                                'accept' = 'RP Acceptor', 'reject' = 'RP Rejector', 'other' = 'Other')) %>%
  rename('Pairing' = 'pairing', 'Type' = 'type',
         'Environment' = 'environment', 'Count' = 'count') %>%
  mutate_at(c('Count', 'Percent'), as.numeric) %>%
  mutate(Percent = format(round(Percent, 2), nsmall = 2)) %>%
  mutate_all(as.character())

dt %>% write_rds('tables/ajt_acceptor.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/ajt_acceptor.csv')
dt <- read_rds('tables/ajt_acceptor.rds')

# ---------------------------------------------------------------------------- #
# model: acceptor ~ environment * pairing ####
# ---------------------------------------------------------------------------- #

md <- ajt_ind %>% 
  mutate(acceptor = case_when(type == 'accept' ~ TRUE, TRUE ~ FALSE)) %>%
  select(study, pairing, participant, environment, acceptor) %>%
  mutate_at(c('study', 'pairing', 'participant', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

# apply deviation coding ...

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$pairing) <- contr.treatment(5) - matrix(rep(1/5, 20), ncol = 4)

# view contrasts ...

contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- glmer(acceptor ~ environment * pairing + 
                  (1 | participant),
                data = filter(md, study == 'ORC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 4 sec
# singular fit

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/orc_ajt_acceptor_model.rds')
model <- read_rds('models/orc_ajt_acceptor_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_acceptor_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_acceptor_model.csv')
dt <- read_rds('tables/orc_ajt_acceptor_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                       |Est.  |$SE$ |$z$   |$p$    |
# |:-------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                     |−1.86 |0.20 |−9.09 |< .001 | ***
# |Environment[Long]               |0.66  |0.26 |2.57  |0.010  | *
# |Environment[Island]             |1.06  |0.24 |4.37  |< .001 | ***
# |Pairing[KE]                     |2.83  |0.54 |5.23  |< .001 | ***
# |Pairing[ME]                     |3.49  |0.54 |6.44  |< .001 | ***
# |Pairing[KK]                     |2.32  |0.54 |4.26  |< .001 | ***
# |Pairing[MM]                     |2.54  |0.54 |4.74  |< .001 | ***
# |Environment[Long]:Pairing[KE]   |0.45  |0.96 |0.47  |0.639  |
# |Environment[Island]:Pairing[KE] |−0.42 |0.85 |−0.50 |0.616  |
# |Environment[Long]:Pairing[ME]   |0.00  |0.94 |0.00  |> .999 |
# |Environment[Island]:Pairing[ME] |−2.05 |0.83 |−2.47 |0.014  | *
# |Environment[Long]:Pairing[KK]   |1.34  |1.00 |1.35  |0.178  |
# |Environment[Island]:Pairing[KK] |0.09  |0.88 |0.10  |0.917  |
# |Environment[Long]:Pairing[MM]   |1.51  |0.97 |1.56  |0.119  |
# |Environment[Island]:Pairing[MM] |−1.17 |0.86 |−1.36 |0.174  |

# joint tests ...

joint_tests(model)

# model term          df1 df2 F.ratio p.value
# environment           2 Inf   9.561  0.0001 ***
# pairing               4 Inf  12.245  <.0001 ***
# environment:pairing   8 Inf   3.002  0.0023 **

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ environment * pairing) %>%
  contrast('revpairwise', simple = 'environment', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast       pairing  estimate    SE  df z.ratio p.value
# long - short   EE       6.40e-06 0.837 Inf   0.000  1.0000
# island - short EE       1.77e+00 0.710 Inf   2.499  0.1551
# island - long  EE       1.77e+00 0.710 Inf   2.499  0.1551
# long - short   KE       4.52e-01 0.477 Inf   0.947  0.9909
# island - short KE       1.35e+00 0.470 Inf   2.872  0.0553 .
# island - long  KE       8.97e-01 0.452 Inf   1.983  0.4599
# long - short   ME      -1.80e-06 0.424 Inf   0.000  1.0000
# island - short ME      -2.73e-01 0.427 Inf  -0.639  0.9997
# island - long  ME      -2.73e-01 0.427 Inf  -0.639  0.9997
# long - short   KK       1.34e+00 0.540 Inf   2.488  0.1598
# island - short KK       1.86e+00 0.535 Inf   3.485  0.0071 *
# island - long  KK       5.22e-01 0.459 Inf   1.137  0.9657
# long - short   MM       1.51e+00 0.484 Inf   3.111  0.0261 *
# island - short MM       6.02e-01 0.495 Inf   1.217  0.9459
# island - long  MM      -9.05e-01 0.455 Inf  -1.988  0.4562

# table ...

dt <- pairwise %>%
  select(pairing, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_acceptor_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_acceptor_posthoc.csv')
dt <- read_rds('tables/orc_ajt_acceptor_posthoc.rds')

dt %>% rename("Est." = 3) %>% kable(format = "markdown")

# |Pairing |Contrast       |Est.  |$SE$ |$z$   |$p$    |
# |:-------|:--------------|:-----|:----|:-----|:------|
# |EE      |Long − Short   |0.00  |0.84 |0.00  |> .999 |
# |EE      |Island − Short |1.77  |0.71 |2.50  |0.155  |
# |EE      |Island − Long  |1.77  |0.71 |2.50  |0.155  |
# |KE      |Long − Short   |0.45  |0.48 |0.95  |0.991  |
# |KE      |Island − Short |1.35  |0.47 |2.87  |0.055  | .
# |KE      |Island − Long  |0.90  |0.45 |1.98  |0.460  |
# |ME      |Long − Short   |0.00  |0.42 |0.00  |> .999 |
# |ME      |Island − Short |−0.27 |0.43 |−0.64 |> .999 |
# |ME      |Island − Long  |−0.27 |0.43 |−0.64 |> .999 |
# |KK      |Long − Short   |1.34  |0.54 |2.49  |0.160  |
# |KK      |Island − Short |1.86  |0.54 |3.48  |0.007  | **
# |KK      |Island − Long  |0.52  |0.46 |1.14  |0.966  |
# |MM      |Long − Short   |1.51  |0.48 |3.11  |0.026  | *
# |MM      |Island − Short |0.60  |0.49 |1.22  |0.946  |
# |MM      |Island − Long  |−0.90 |0.46 |−1.99 |0.456  |

# plot ...

plot(emmeans(model, ~ environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- glmer(acceptor ~ environment * pairing + 
                  (1 | participant),
                data = filter(md, study == 'SRC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 5 sec
# no warnings

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/src_ajt_acceptor_model.rds')
model <- read_rds('models/src_ajt_acceptor_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_ajt_acceptor_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_acceptor_model.csv')
dt <- read_rds('tables/src_ajt_acceptor_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                       |Est.  |$SE$ |$z$   |$p$    |
# |:-------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                     |−1.26 |0.15 |−8.39 |< .001 | ***
# |Environment[Long]               |1.20  |0.35 |3.47  |< .001 | ***
# |Environment[Island]             |1.74  |0.29 |5.98  |< .001 | ***
# |Pairing[KE]                     |2.59  |0.57 |4.57  |< .001 | ***
# |Pairing[ME]                     |4.00  |0.57 |7.06  |< .001 | ***
# |Pairing[KK]                     |2.45  |0.58 |4.20  |< .001 | ***
# |Pairing[MM]                     |3.22  |0.56 |5.72  |< .001 | ***
# |Environment[Long]:Pairing[KE]   |0.82  |1.51 |0.54  |0.587  |
# |Environment[Island]:Pairing[KE] |−0.82 |1.18 |−0.69 |0.488  |
# |Environment[Long]:Pairing[ME]   |0.82  |1.49 |0.55  |0.583  |
# |Environment[Island]:Pairing[ME] |−1.99 |1.15 |−1.74 |0.082  | .
# |Environment[Long]:Pairing[KK]   |2.98  |1.56 |1.90  |0.057  | .
# |Environment[Island]:Pairing[KK] |0.37  |1.24 |0.29  |0.768  |
# |Environment[Long]:Pairing[MM]   |1.29  |1.49 |0.87  |0.385  |
# |Environment[Island]:Pairing[MM] |−1.23 |1.15 |−1.06 |0.287  |

# joint tests ...

joint_tests(model)

# model term          df1 df2 F.ratio p.value
# environment           2 Inf  18.120  <.0001 ***
# pairing               4 Inf  15.707  <.0001 ***
# environment:pairing   8 Inf   2.751  0.0049 **

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ environment * pairing) %>%
  contrast('revpairwise', simple = 'environment', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast       pairing estimate    SE  df z.ratio p.value
# long - short   EE        0.0161 1.437 Inf   0.011  1.0000
# island - short EE        2.4712 1.086 Inf   2.276  0.2591
# island - long  EE        2.4551 1.085 Inf   2.262  0.2672
# long - short   KE        0.8369 0.467 Inf   1.793  0.6050
# island - short KE        1.6560 0.457 Inf   3.620  0.0042 **
# island - long  KE        0.8190 0.410 Inf   1.995  0.4465
# long - short   ME        0.8322 0.377 Inf   2.206  0.3003
# island - short ME        0.4782 0.371 Inf   1.288  0.9211
# island - long  ME       -0.3541 0.377 Inf  -0.939  0.9913
# long - short   KK        2.9921 0.614 Inf   4.873  <.0001 ***
# island - short KK        2.8380 0.613 Inf   4.627  0.0001 ***
# island - long  KK       -0.1540 0.393 Inf  -0.392  1.0000
# long - short   MM        1.3107 0.395 Inf   3.320  0.0127 *
# island - short MM        1.2435 0.394 Inf   3.152  0.0225 *
# island - long  MM       -0.0672 0.367 Inf  -0.183  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_acceptor_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_acceptor_posthoc.csv')
dt <- read_rds('tables/src_ajt_acceptor_posthoc.rds')

dt %>% rename("Est." = 3) %>% kable(format = "markdown")

# |Pairing |Contrast       |Est.  |$SE$ |$z$   |$p$    |
# |:-------|:--------------|:-----|:----|:-----|:------|
# |EE      |Long − Short   |0.02  |1.44 |0.01  |> .999 |
# |EE      |Island − Short |2.47  |1.09 |2.28  |0.259  |
# |EE      |Island − Long  |2.46  |1.09 |2.26  |0.267  |
# |KE      |Long − Short   |0.84  |0.47 |1.79  |0.605  |
# |KE      |Island − Short |1.66  |0.46 |3.62  |0.004  | **
# |KE      |Island − Long  |0.82  |0.41 |2.00  |0.447  |
# |ME      |Long − Short   |0.83  |0.38 |2.21  |0.300  |
# |ME      |Island − Short |0.48  |0.37 |1.29  |0.921  |
# |ME      |Island − Long  |−0.35 |0.38 |−0.94 |0.991  |
# |KK      |Long − Short   |2.99  |0.61 |4.87  |< .001 | ***
# |KK      |Island − Short |2.84  |0.61 |4.63  |< .001 | ***
# |KK      |Island − Long  |−0.15 |0.39 |−0.39 |> .999 |
# |MM      |Long − Short   |1.31  |0.39 |3.32  |0.013  | *
# |MM      |Island − Short |1.24  |0.39 |3.15  |0.022  | *
# |MM      |Island − Long  |−0.07 |0.37 |−0.18 |> .999 |

# plot ...

plot(emmeans(model, ~ environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: proficiency * category * environment (bar) ####
# ---------------------------------------------------------------------------- #

temp <- ajt_ind %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(is.na(proficiency) == FALSE, is.na(type) == FALSE) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study) %>%
  mutate(prof = scale(proficiency)) %>%
  ungroup() %>%
  select(study, group, group_label, participant, environment, proficiency, prof, type) %>%
  mutate(type = fct_relevel(type, "accept", "undefined", "reject"))

plot <- temp %>%
  group_by(study, group, environment, type) %>%
  mutate(n = n()) %>%
  ungroup()

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list( 
  aes(x=environment, y=prof, col=type, fill=type),
  geom_hline(yintercept = 0),
  #geom_violin(width = .1, position = position_dodge(width = .5)),
  geom_boxplot(outlier.size = .5,
               width = .25,
               fill = 'white',
               position = position_dodge(width = .5, preserve = 'single')),
  theme_classic(),
  coord_cartesian(ylim = c(-2.85, 2.85)),
  scale_x_discrete(name = 'Environment',
                   limits = c('short', 'long', 'island'),
                   labels=c('Short', 'Long', 'Island')),
  scale_y_continuous(name='Proficiency (*z*-score)'),
  scale_fill_manual(name='Participant Type',
                    values=c('#b0e2de', '#b5dff5', '#ffccb9'),
                    limits = c('accept', 'undefined', 'reject'),
                    labels=c('RP Acceptor','Inconsistent', 'RP Rejector')),
  scale_colour_manual(name='Participant Type',
                      values=c('#40B0A6', '#47b0e5', 'coral'),
                      limits = c('accept', 'undefined', 'reject'),
                      labels=c('RP Acceptor', 'Inconsistent', 'RP Rejector')),
  theme(text = element_text(size = 12),
        axis.title.y = element_markdown(),
        strip.text = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = 'right'),
  facet_grid2(~group_label, axes = 'all', remove_labels = FALSE) 
)

p1 + s 
ggsave("plots/orc_ajt_acceptor_environment_proficiency.png", width=6.5, height=2.5, dpi=600)

p2 + s 
ggsave("plots/src_ajt_acceptor_environment_proficiency.png", width=6.5, height=2.5, dpi=600)

# ---------------------------------------------------------------------------- #
# model: acceptor ~ environment * group * proficiency ####
# ---------------------------------------------------------------------------- #

md <- plot %>%
  mutate(acceptor = case_when(type == 'accept' ~ TRUE, TRUE ~ FALSE)) %>%
  select(study, group, participant, environment, proficiency = prof, acceptor) %>%
  mutate_at(c('study', 'group', 'participant', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

# apply deviation coding ...

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)

# view contrasts ...

contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- glmer(acceptor ~ environment * group * proficiency + 
                  (1 | participant),
                data = filter(md, study == 'ORC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 1 sec
# no warnings

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/orc_ajt_acceptor_proficiency_model.rds')
model <- read_rds('models/orc_ajt_acceptor_proficiency_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_acceptor_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_acceptor_proficiency_model.csv')
dt <- read_rds('tables/orc_ajt_acceptor_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                  |Est.  |$SE$ |$z$   |$p$    |
# |:------------------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                                |−1.47 |0.40 |−3.73 |< .001 | ***
# |Environment[Long]                          |0.40  |0.39 |1.01  |0.313  |
# |Environment[Island]                        |0.83  |0.40 |2.05  |0.040  | *
# |Group[KLE]                                 |0.94  |0.68 |1.39  |0.166  |
# |Proficiency                                |−0.75 |0.36 |−2.10 |0.036  | *
# |Environment[Long]:Group[KLE]               |−0.83 |0.79 |−1.05 |0.294  |
# |Environment[Island]:Group[KLE]             |−2.51 |0.83 |−3.03 |0.002  | **
# |Environment[Long]:Proficiency              |0.21  |0.39 |0.54  |0.590  |
# |Environment[Island]:Proficiency            |0.15  |0.40 |0.38  |0.706  |
# |Group[KLE]:Proficiency                     |0.37  |0.67 |0.55  |0.583  |
# |Environment[Long]:Group[KLE]:Proficiency   |−0.84 |0.79 |−1.06 |0.287  |
# |Environment[Island]:Group[KLE]:Proficiency |−1.12 |0.81 |−1.39 |0.164  |

# joint tests ...

joint_tests(model)

# model term                    df1 df2 F.ratio p.value
# environment                     2 Inf   2.117  0.1204
# group                           1 Inf   1.918  0.1660
# proficiency                     1 Inf   4.392  0.0361 *
# environment:group               2 Inf   4.871  0.0077 *
# environment:proficiency         2 Inf   0.150  0.8605
# group:proficiency               1 Inf   0.301  0.5832
# environment:group:proficiency   2 Inf   1.033  0.3559

# post-hoc tests ...

pairwise <- test(emtrends(model, ~ environment * group, var = "proficiency"))
pairwise

# environment group proficiency.trend    SE  df z.ratio p.value
# short       KLE              -1.383 0.695 Inf  -1.989  0.0467 *
# long        KLE              -0.752 0.632 Inf  -1.191  0.2338
# island      KLE              -0.671 0.602 Inf  -1.114  0.2651
# short       MLE              -0.360 0.513 Inf  -0.701  0.4830
# long        MLE              -0.567 0.518 Inf  -1.096  0.2731
# island      MLE              -0.771 0.532 Inf  -1.450  0.1472

# table ...

dt <- pairwise %>%
  select(group, environment, Trend = proficiency.trend, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_acceptor_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_acceptor_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_ajt_acceptor_proficiency_posthoc.rds')

dt %>% rename("Est." = 3) %>% kable(format = "markdown")

# |Group |Environment |Est.  |$SE$ |$z$   |$p$   |
# |:-----|:-----------|:-----|:----|:-----|:-----|
# |KLE   |Short       |−1.38 |0.70 |−1.99 |0.047 | *
# |KLE   |Long        |−0.75 |0.63 |−1.19 |0.234 |
# |KLE   |Island      |−0.67 |0.60 |−1.11 |0.265 |
# |MLE   |Short       |−0.36 |0.51 |−0.70 |0.483 |
# |MLE   |Long        |−0.57 |0.52 |−1.10 |0.273 |
# |MLE   |Island      |−0.77 |0.53 |−1.45 |0.147 |

# plot ...

plot(emmeans(model, ~ environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- glmer(acceptor ~ environment * group * proficiency + 
                  (1 | participant),
                data = filter(md, study == 'SRC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 1 sec
# no warnings

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/src_ajt_acceptor_proficiency_model.rds')
model <- read_rds('models/src_ajt_acceptor_proficiency_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_ajt_acceptor_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_acceptor_proficiency_model.csv')
dt <- read_rds('tables/src_ajt_acceptor_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                  |Est.  |$SE$ |$z$   |$p$    |
# |:------------------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                                |−0.48 |0.19 |−2.50 |0.012  | *
# |Environment[Long]                          |1.07  |0.34 |3.09  |0.002  | **
# |Environment[Island]                        |1.25  |0.34 |3.64  |< .001 | ***
# |Group[KLE]                                 |1.63  |0.40 |4.07  |< .001 | ***
# |Proficiency                                |−0.39 |0.20 |−1.99 |0.046  | *
# |Environment[Long]:Group[KLE]               |0.11  |0.67 |0.16  |0.873  |
# |Environment[Island]:Group[KLE]             |−1.15 |0.67 |−1.72 |0.086  | .
# |Environment[Long]:Proficiency              |0.54  |0.35 |1.56  |0.120  |
# |Environment[Island]:Proficiency            |0.92  |0.35 |2.63  |0.008  | **
# |Group[KLE]:Proficiency                     |0.03  |0.38 |0.07  |0.945  |
# |Environment[Long]:Group[KLE]:Proficiency   |0.31  |0.69 |0.45  |0.655  |
# |Environment[Island]:Group[KLE]:Proficiency |−0.65 |0.69 |−0.95 |0.343  |

# joint tests ...

joint_tests(model)

# model term                    df1 df2 F.ratio p.value
# environment                     2 Inf   7.291  0.0007 ***
# group                           1 Inf  16.551  <.0001 ***
# proficiency                     1 Inf   3.966  0.0464 *
# environment:group               2 Inf   2.405  0.0903 .
# environment:proficiency         2 Inf   3.471  0.0311 *
# group:proficiency               1 Inf   0.005  0.9449
# environment:group:proficiency   2 Inf   1.175  0.3089

# post-hoc tests ...

pairwise <- test(emtrends(model, ~ environment * group, var = "proficiency"))
pairwise

# environment group proficiency.trend    SE  df z.ratio p.value
# short       KLE              -0.948 0.471 Inf  -2.012  0.0442 *
# long        KLE              -0.561 0.405 Inf  -1.385  0.1661
# island      KLE               0.302 0.384 Inf   0.786  0.4319
# short       MLE              -0.807 0.369 Inf  -2.186  0.0288 *
# long        MLE              -0.110 0.349 Inf  -0.316  0.7522
# island      MLE              -0.210 0.343 Inf  -0.613  0.5402

# table ...

dt <- pairwise %>%
  select(group, environment, Trend = proficiency.trend, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_acceptor_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_acceptor_proficiency_posthoc.csv')
dt <- read_rds('tables/src_ajt_acceptor_proficiency_posthoc.rds')

dt %>% rename("Est." = 3) %>% kable(format = "markdown")

# |Group |Environment |Est.  |$SE$ |$z$   |$p$   |
# |:-----|:-----------|:-----|:----|:-----|:-----|
# |KLE   |Short       |−0.95 |0.47 |−2.01 |0.044 | *
# |KLE   |Long        |−0.56 |0.40 |−1.38 |0.166 |
# |KLE   |Island      |0.30  |0.38 |0.79  |0.432 |
# |MLE   |Short       |−0.81 |0.37 |−2.19 |0.029 | *
# |MLE   |Long        |−0.11 |0.35 |−0.32 |0.752 |
# |MLE   |Island      |−0.21 |0.34 |−0.61 |0.540 |

# plot ...

plot(emmeans(model, ~ environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: proficiency * category * environment (line) ####
# ---------------------------------------------------------------------------- #

# summarize ...

plot2 <- plot %>%
  group_by(study, group, group_label, type, environment) %>%
  summarise(mean = mean(prof, na.rm = T),
            lcl = my.lcl(prof),
            ucl = my.ucl(prof)) %>%
  ungroup() %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island")) %>%
  arrange(study, group, type, environment)

# plot ...

p1 <- ggplot(data=filter(plot2, study == 'ORC'))
p2 <- ggplot(data=filter(plot2, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=type, col=type, 
      shape=type),
  geom_hline(yintercept = 0),
  geom_line(aes(linetype = type), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                width = 0.2, 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name='Environment'),
  scale_y_continuous(name='Proficiency (*z*-score)'),
  coord_cartesian(ylim = c(-2, 2)),
  scale_colour_manual(name='Category', 
                      values=c('#40B0A6', 'coral', '#47b0e5'), 
                      limits = c('accept', 'reject', 'undefined'),
                      labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  scale_linetype_manual(name='Category', 
                        values=c('solid', '41', '7222'), 
                        limits = c('accept', 'reject', 'undefined'),
                        labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        axis.title.y = element_markdown(),
        legend.position = 'right',
        legend.margin = margin(t = -.4, unit = 'cm'),
        legend.key.width = unit(.75, "cm"),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~group_label, design = 'AB', scales = 'free')
)

p1 + s 
ggsave("plots/orc_ajt_acceptor_environment_proficiency2.png", width=6.5, height=2.25, dpi=600)

p2 + s 
ggsave("plots/src_ajt_acceptor_environment_proficiency2.png", width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot2 %>%
  select(study, group, environment, mean, lcl, ucl) %>%
  arrange(study, group, environment) %>%
  my_table() 

dt %>% write_rds('tables/ajt_acceptor_environment_proficiency2.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/ajt_acceptor_environment_proficiency2.csv')
dt <- read_rds('tables/ajt_acceptor_environment_proficiency2.rds')

# ---------------------------------------------------------------------------- #
# plot: RP acceptors vs. non-acceptors vs. inconsistent raters (bar) ####
# ---------------------------------------------------------------------------- #

# plot by group ...

ajt_ind2 <- ajt_ind %>%
  group_by(study, group, language, pairing, participant, proficiency) %>%
  summarise(max = max(acceptance)) %>%
  ungroup() %>%
  mutate(subgroup = case_when(max >= .8 ~ "RP Acceptor",
                              max <= .2 ~ "RP Rejector",
                              TRUE ~ "Inconsistent"))

plot <- ajt_ind2 %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, subgroup) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(study, pairing) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(prop = count / total) %>%
  mutate(subgroup = fct_relevel(subgroup, 
                                "RP Acceptor", 
                                "Inconsistent",
                                "RP Rejector" 
  )) %>%
  mutate(pairing_num = case_when(pairing == "EE" ~ 1,
                                 pairing == "KE" ~ 2,
                                 pairing == "ME" ~ 3,
                                 pairing == "KK" ~ 4,
                                 pairing == "MM" ~ 5)) %>%
  mutate(pairing_label = fct_reorder(pairing_label, pairing_num))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=pairing_label, y=(count/total)*100, group=subgroup, fill=subgroup, col=subgroup, 
      label=paste0(round(count/total * 100), '% (*n* = ', count, ')')),
  geom_bar(stat = "identity", width = .85, 
           col = 'black', linewidth = .25
           ),
  geom_richtext(col = 'black', fill = NA, label.color = NA, size = 3,
                position = position_stack(vjust = 0.5)),
  theme_classic(),
  scale_x_discrete(name='Pairing'),
  scale_y_continuous(name='Participants (%)', 
                     #limits=c(0, 101), 
                     #breaks=c(0, 20, 40, 60, 80, 100)
  ),
  scale_fill_manual(name="Type:", 
                    breaks = c('RP Acceptor', 'RP Rejector', 'Inconsistent'),
                    values=c('#b0e2de', '#ffccb9', '#b5dff5')),
  # scale_color_manual(name="Type:",
  #                    values=c('#40B0A6', '#47b0e5', 'coral')),
  theme(text = element_text(size = 12),
        #axis.title.x = element_blank(),
        #legend.title = element_blank(),
        axis.text.x = element_markdown(),
        legend.box.margin = margin(-10, 0, 0, 0),
        legend.position = "bottom", 
        legend.margin=margin(1, 1, 1, 1))
)

p1 + s 
ggsave('plots/orc_ajt_subgroup_prop.png', width=6.5, height=3.5, dpi=600)

p2 + s
ggsave('plots/src_ajt_subgroup_prop.png', width=6.5, height=3.5, dpi=600)

# plot by group and environment ...

plot <- ajt_ind %>%
  mutate(subgroup = case_when(type == 'accept' ~ 'RP Acceptor',
                              type == 'reject' ~ 'RP Rejector',
                              type == 'undefined' ~ 'Inconsistent')) %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, subgroup, environment) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(study, pairing, environment) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(prop = count / total) %>%
  mutate(subgroup = fct_relevel(subgroup, 
                                "RP Acceptor", 
                                "Inconsistent",
                                "RP Rejector"),
         environment = fct_relevel(environment,
                                   "short",
                                   "long",
                                   "island"))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=(count/total)*100, group=subgroup, fill=subgroup, col=subgroup, 
      label=paste0(round(count/total * 100))),
  geom_bar(stat = "identity", width = .9, col = 'black', linewidth = .2),
  geom_richtext(col = 'black', fill = NA, label.color = NA, size = 3,
                position = position_stack(vjust = 0.5)),
  theme_classic(),
  scale_x_discrete(name='Environment',
                   breaks = c('short', 'long', 'island'),
                   labels = c('Short', 'Long', 'Island')),
  scale_y_continuous(name='Participants (%)'
                     #breaks=c(0, 20, 40, 60, 80, 100)
  ),
  scale_fill_manual(name="Type:", 
                    breaks = c('RP Acceptor', 'Inconsistent', 'RP Rejector'),
                    values=c('#b0e2de', '#b5dff5', '#ffccb9')),
  scale_color_manual(name="Type:",
                     values=c('#40B0A6', '#47b0e5', 'coral'),
                     breaks = c('RP Acceptor', 'Inconsistent', 'RP Rejector')),
  theme(text = element_text(size = 12),
        #axis.title.x = element_blank(),
        #legend.title = element_blank(),
        strip.text = element_markdown(),
        panel.spacing = unit(.1, "lines"),
        legend.box.margin = margin(-10, 0, 0, 0),
        legend.position = "bottom", 
        legend.margin=margin(1, 1, 1, 1)),
  facet_manual(~pairing_label, design = 'ABDCE')
)

p1 + s
ggsave('plots/orc_ajt_subgroup_prop_environment.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_ajt_subgroup_prop_environment.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: RP acceptors by type ####
# ---------------------------------------------------------------------------- #

# list RP acceptors by environment in L1 and L2

temp <- ajt_ind %>%
  filter(pairing != 'EE') %>%
  mutate(acceptor = case_when(acceptance >= .8 ~ TRUE, TRUE ~ FALSE),
         lang = case_when(pairing %in% c('KE', 'ME') ~ 'L2',
                          pairing %in% c('KK', 'MM') ~ 'L1')) %>%
  group_by(study, group, lang, participant, environment, acceptor) %>%
  summarise() %>%
  ungroup() %>%
  pivot_wider(names_from = lang, values_from = acceptor) %>%
  mutate(transfer = case_when(L1 > 0 & L2 > 0 ~ TRUE, TRUE ~ FALSE))

list1 <- temp %>%
  group_by(study, group, participant) %>%
  summarise(acceptor = mean(L2),
            transfer = mean(transfer)) %>%
  ungroup() %>%
  mutate(acceptor = case_when(acceptor > 0 ~ TRUE, TRUE ~ FALSE),
         transfer = case_when(transfer > 0 ~ TRUE, TRUE ~ FALSE))

# list yes bias participants

temp1 <- ajt_clean %>% 
  filter(pairing %in% c('KE', 'ME')) %>%
  group_by(study, group, participant) %>%
  summarise(rate = mean(acceptance)) %>%
  ungroup() %>%
  mutate(yesbias1 = case_when(rate >= .8 ~ TRUE, TRUE ~ FALSE)) %>%
  select(study, group, participant, yesbias1)

temp2 <- ajt_clean %>% 
  filter(pairing %in% c('KE', 'ME'), condition == 'ungrammatical') %>%
  group_by(study, group, participant) %>%
  summarise(rate = mean(acceptance)) %>%
  ungroup() %>%
  mutate(yesbias2 = case_when(rate >= .8 ~ TRUE, TRUE ~ FALSE)) %>%
  select(study, group, participant, yesbias2)

list2 <- temp1 %>%
  left_join(temp2, by = c('study', 'group', 'participant')) %>%
  mutate(yesbias = case_when(yesbias1 == TRUE & yesbias2 == TRUE ~ TRUE, TRUE ~ FALSE)) %>%
  select(-yesbias1, -yesbias2)

# combine and plot ...

temp <- list1 %>%
  left_join(list2, by = c('study', 'group', 'participant')) %>%
  mutate(neither = case_when(transfer == FALSE & yesbias == FALSE ~ TRUE, 
                             TRUE ~ FALSE))

temp2 <- temp %>%
  mutate(type = case_when(yesbias == TRUE ~ "yesbias",
                          transfer == TRUE ~ "transfer",
                          neither == TRUE ~ "neither"))

plot <- temp2 %>%
  filter(acceptor == TRUE) %>%
  group_by(study, group) %>%
  summarise(total = length(unique(participant)),
            transfer = length(which(type == "transfer")),
            yesbias = length(which(type == "yesbias")),
            neither = length(which(type == "neither"))) %>%
  pivot_longer(cols = c(transfer:neither), names_to = 'type', values_to = 'n') %>%
  ungroup() %>%
  mutate(prop = n / total) %>%
  mutate(type = fct_relevel(type, 'transfer', 'yesbias', 'neither'),
         group_label = paste0(group, '<br>(*n* = ', total, ')')) %>%
  mutate(type = fct_relevel(type, "yesbias", "transfer", "neither"))
  
p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=type, y=prop * 100, group=type,
      label = paste0(round(prop * 100), '% (*n* = ', n, ')')),
  geom_richtext(col = 'black', hjust=.5, vjust = -.1, size=3,
                fill = NA, label.color = NA),
  geom_bar(stat = "identity", width = .85, fill='#b0e2de', col='black',
           linewidth = .25),
  theme_classic(),
  coord_cartesian(ylim = c(0, 100)),
  scale_x_discrete(name = 'RP Acceptor Type',
                   breaks = c('yesbias', 'transfer', 'neither'),
                   labels=c('Bias', 'Transfer', 'Other')),
  scale_y_continuous(name='Percentage', breaks = c(0, 50, 100)),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        axis.text.x = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = 'hide'
  ),
  facet_grid2(~group_label, axes = 'all', remove_labels = FALSE) 
)

p1 + s 
ggsave("plots/orc_ajt_acceptor_type.png", width=6.5, height=2, dpi=600)

p2 + s 
ggsave("plots/src_ajt_acceptor_type.png", width=6.5, height=2, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: proficiency (acceptor vs. rejector vs. inconsistent) ####
# ---------------------------------------------------------------------------- #

temp <- ajt_ind %>%
  group_by(study, group, language, pairing, participant, proficiency) %>%
  summarise(max = max(acceptance)) %>%
  ungroup() %>%
  mutate(type = case_when(max >= .8 ~ "acceptor",
                          max <= .2 ~ "rejector",
                          TRUE ~ "inconsistent"))

check <- temp %>%
  filter(study == 'SRC') %>%
  filter(language == 'English') %>%
  group_by(study, pairing, type) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(study, pairing) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  arrange(study, type, pairing)

temp2 <- temp %>% 
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(is.na(proficiency) == FALSE, is.na(type) == FALSE) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study) %>%
  mutate(prof = scale(proficiency)) %>%
  ungroup() %>%
  select(study, group, group_label, participant, proficiency, prof, type) %>%
  mutate(type = fct_relevel(type, "acceptor", "rejector", "inconsistent"))

plot <- temp2 %>%
  group_by(study, group, type) %>%
  mutate(n = n()) %>%
  ungroup()
  
p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=type, y=prof, group=type, col=type, fill=type, shape=type, label=n),
  geom_hline(yintercept = 0),
  geom_violin(),
  geom_richtext(aes(label = paste0('*n* = ', n), y = -2.75), 
                col = 'black', fill = NA, label.color = NA, size = 3.3),
  geom_boxplot(width = .1, fill = 'white'),
  theme_classic(),
  coord_cartesian(ylim = c(-2.85, 2.85)),
  scale_x_discrete(name = 'Participant Type',
                   limits = c('acceptor', 'rejector', 'inconsistent'),
                   labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  scale_y_continuous(name='Proficiency (*z*-score)'),
  scale_fill_manual(name='Participant Type', 
                      values=c('#b0e2de', '#ffccb9', '#b5dff5'),
                      limits = c('acceptor', 'rejector', 'inconsistent'),
                      labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  scale_colour_manual(name='Participant Type', 
                      values=c('#40B0A6', 'coral', '#47b0e5'),
                      limits = c('acceptor', 'rejector', 'inconsistent'),
                      labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  theme(text = element_text(size = 12),
        axis.title.y = element_markdown(),
        strip.text = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = 'hide'
  ),
  facet_grid2(~group_label, axes = 'all', remove_labels = FALSE) 
)

p1 + s 
ggsave("plots/orc_ajt_acceptor_type_proficiency.png", width=6.5, height=2.5, dpi=600)

p2 + s 
ggsave("plots/src_ajt_acceptor_type_proficiency.png", width=6.5, height=2.5, dpi=600)

# ---------------------------------------------------------------------------- #
# model: proficiency (acceptor vs. rejector) ####
# ---------------------------------------------------------------------------- #

md <- temp2 %>% 
  select(-group_label) %>%
  mutate_at(c('study', 'group', 'participant', 'type'), factor) %>%
  mutate(type = fct_relevel(type, c('acceptor', 'rejector'))) %>%
  ungroup()

summary(md)

# apply deviation coding ...

contrasts(md$type) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

# view contrasts ...

contrasts(md$type)
contrasts(md$group)

# fit orc model ...

model <- lm(prof ~ group * type, 
               data = filter(md, study == 'ORC'))
summary(model)

# Parameter    Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   0.07089    0.09014   0.786 0.432963    
# group2       -0.03724    0.18028  -0.207 0.836668    
# type2         0.65292    0.19121   3.415 0.000839 ***
# type3         0.22810    0.22147   1.030 0.304859    
# group2:type2 -0.07408    0.38242  -0.194 0.846689    
# group2:type3 -0.25082    0.44295  -0.566 0.572143    

# save ...

model %>% write_rds("models/orc_ajt_acceptor_type_proficiency_model.rds")
model <- read_rds("models/orc_ajt_acceptor_type_proficiency_model.rds")

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('group2' = 'group3',
              'type2' = 'Type[Rejector]',
              'type3' = 'Type[Inconsistent]')) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_acceptor_type_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_acceptor_type_proficiency_model.csv')
dt <- read_rds('tables/orc_ajt_acceptor_type_proficiency_model.rds')

# fit src model ...

model <- lm(prof ~ type, 
            data = filter(md, study == 'SRC'))
summary(model)

# Parameter    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.03852    0.10795   0.357    0.722
# type2        0.23533    0.24715   0.952    0.343
# type3       -0.03368    0.23348  -0.144    0.886

# save ...

model %>% write_rds("models/src_ajt_acceptor_type_proficiency_model.rds")
model <- read_rds("models/src_ajt_acceptor_type_proficiency_model.rds")

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('type2' = 'Type[Rejector]',
              'type3' = 'Type[Inconsistent]')) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_acceptor_type_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_acceptor_type_proficiency_model.csv')
dt <- read_rds('tables/src_ajt_acceptor_type_proficiency_model.rds')

# ---------------------------------------------------------------------------- #
# model: proficiency vs. consistent RP acceptance for L2ers (environment) ####
# ---------------------------------------------------------------------------- #

# plot ...

plot <- ajt_ind %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, participant, proficiency, environment, type) %>%
  summarise() %>%
  ungroup() %>%
  mutate(environment = fct_rev(environment))%>%
  mutate(environment = fct_recode(environment, 
                                  'Short' = 'short', 
                                  'Long' = 'long', 
                                  'Island' = 'island')) %>%
  filter(is.na(proficiency) == FALSE)

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=proficiency, col = type, fill = type),
  geom_hline(aes(yintercept = 25)),
  geom_violin(col = NA, position = position_dodge(width = .6)),
  geom_boxplot(width = .1, fill = 'white', position = position_dodge(width = .6, preserve = 'single')),
  theme_classic(),
  xlab('Environment'),
  scale_y_continuous(name="Proficiency",
                     breaks = c(0, 25, 50),
                     limits=c(0, 50)),
  scale_colour_manual(name='Type', 
                      values=c('#40B0A6', 'coral', '#ffb000'), 
                      limits = c('accept', 'reject', 'undefined'),
                      labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  scale_fill_manual(name='Type', 
                    values=c('#b0e2de', '#ffccb9', '#ffdf99'),
                    limits = c('accept', 'reject', 'undefined'),
                    labels=c('RP Acceptor', 'RP Rejector', 'Inconsistent')),
  #scale_fill_manual(name = 'Group', 
  #values = c('#44AA99', '#DDCC77', '#CC6677')),
  theme(text = element_text(size = 12),
        legend.position = 'right'),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE
  )
)

p1 + s
ggsave("plots/orc_ajt_acceptor_proficiency.png", width=6.5, height=2, dpi=600)
p2 + s
ggsave("plots/src_ajt_acceptor_proficiency.png", width=6.5, height=2, dpi=600)

md <- plot %>%
  mutate(acceptor = case_when(type == 'accept' ~ TRUE, TRUE ~ FALSE)) %>%
  select(study, group, participant, proficiency, environment, acceptor) %>%
  mutate_at(c('study', 'group', 'participant', 'environment'), factor) %>%
  group_by(study) %>%
  mutate(proficiency = scale(as.numeric(proficiency))) %>%
  ungroup() %>%
  filter(is.na(proficiency) == FALSE)

summary(md)

# apply deviation coding ...

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

# view contrasts ...

contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model <- glmer(acceptor ~ environment * group * proficiency + 
                 (1 | participant), 
               data = filter(md, study == 'ORC'), 
               family = binomial, 
               control = glmerControl(
                 optimizer = "bobyqa", 
                 optCtrl=list(maxfun=1e6)))
summary(model)
toc()
beep()

# Parameter                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                      -1.4732     0.3954  -3.726 0.000195 ***
# environment2                      0.3978     0.3940   1.010 0.312691    
# environment3                      0.8291     0.4041   2.052 0.040200 *  
# group2                            0.9378     0.6771   1.385 0.166042    
# proficiency                      -0.7507     0.3582  -2.096 0.036112 *  
# environment2:group2              -0.8272     0.7885  -1.049 0.294123    
# environment3:group2              -2.5059     0.8274  -3.029 0.002457 ** 
# environment2:proficiency          0.2116     0.3925   0.539 0.589718    
# environment3:proficiency          0.1504     0.3984   0.378 0.705693    
# group2:proficiency                0.3693     0.6729   0.549 0.583187    
# environment2:group2:proficiency  -0.8383     0.7879  -1.064 0.287368    
# environment3:group2:proficiency  -1.1231     0.8061  -1.393 0.163531 

# post-hoc tests ...

emtrends <- emtrends(
  model, ~ group * environment, 
  var = 'proficiency', 
  adjust = 'mvt'
)
test <- test(emtrends)
test

# group environment proficiency.trend    SE  df z.ratio p.value
# KLE   Short                  -1.383 0.695 Inf  -1.989  0.2137
# MLE   Short                  -0.360 0.513 Inf  -0.701  0.9645
# KLE   Long                   -0.752 0.632 Inf  -1.191  0.7252
# MLE   Long                   -0.567 0.518 Inf  -1.096  0.7884
# KLE   Island                 -0.671 0.602 Inf  -1.114  0.7765
# MLE   Island                 -0.771 0.532 Inf  -1.450  0.5383

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model <- glmer(acceptor ~ environment * group * proficiency + 
                 (1 | participant), 
               data = filter(md, study == 'SRC'), 
               family = binomial, 
               control = glmerControl(
                 optimizer = "bobyqa", 
                 optCtrl=list(maxfun=1e6)))
summary(model)
toc()
beep()

# Parameter                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                     -0.48412    0.19348  -2.502 0.012344 *  
# environment2                     1.06545    0.34436   3.094 0.001975 ** 
# environment3                     1.25251    0.34408   3.640 0.000273 ***
# group2                           1.63410    0.40167   4.068 4.74e-05 ***
# proficiency                     -0.38906    0.19536  -1.992 0.046422 *  
# environment2:group2              0.10755    0.67120   0.160 0.872697    
# environment3:group2             -1.14859    0.66943  -1.716 0.086205 .  
# environment2:proficiency         0.54190    0.34838   1.555 0.119829    
# environment3:proficiency         0.92319    0.35041   2.635 0.008425 ** 
# group2:proficiency               0.02662    0.38490   0.069 0.944864    
# environment2:group2:proficiency  0.30978    0.69269   0.447 0.654723    
# environment3:group2:proficiency -0.65236    0.68822  -0.948 0.343183

# post-hoc tests ...

emtrends <- emtrends(
  model, ~ group * environment, 
  var = 'proficiency', 
  adjust = 'mvt'
)
test <- test(emtrends)
test

# group environment proficiency.trend    SE  df z.ratio p.value
# KLE   Short                  -0.948 0.471 Inf  -2.012  0.2309
# MLE   Short                  -0.807 0.369 Inf  -2.186  0.1568
# KLE   Long                   -0.561 0.405 Inf  -1.385  0.6481
# MLE   Long                   -0.110 0.349 Inf  -0.316  0.9997
# KLE   Island                  0.302 0.384 Inf   0.786  0.9623
# MLE   Island                 -0.210 0.343 Inf  -0.613  0.9892

# ---------------------------------------------------------------------------- #
# model: proficiency vs. consistent RP acceptance for L2ers (overall) ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

plot2 <- plot %>%
  mutate(acceptance = case_when(type == 'accept' ~ TRUE, TRUE ~ FALSE)) %>%
  group_by(study, group, participant, proficiency) %>%
  summarise(acceptance = mean(acceptance)) %>%
  ungroup() %>%
  mutate(acceptor = case_when(acceptance > 0 ~ TRUE, TRUE ~ FALSE)) %>%
  select(-acceptance) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  mutate_at(c('study', 'group', 'participant'), factor) %>%
  group_by(study, group, acceptor) %>%
  mutate(acceptor_size = length(unique(participant))) %>%
  ungroup()

p1 <- ggplot(data=filter(plot2, study == 'ORC'))
p2 <- ggplot(data=filter(plot2, study == 'SRC'))

s <- list(
  aes(x=acceptor, y=proficiency, col = acceptor, fill = acceptor),
  geom_hline(aes(yintercept = 25)),
  geom_text(y = 3, aes(label = paste0('n = ', acceptor_size)),
            hjust=.5, size=3),
  geom_violin(col = NA, position = position_dodge(width = .6)),
  geom_boxplot(width = .1, fill = 'white', position = position_dodge(width = .6, preserve = 'single')),
  theme_classic(),
  scale_x_discrete(name = 'Type',
                   limits = c(TRUE, FALSE), 
                   labels = c('RP Acceptor', 'RP Non-Acceptor')),
  scale_y_continuous(name='Proficiency',
                     breaks = c(0, 25, 50),
                     limits=c(0, 50)),
  scale_colour_manual(name='Type', 
                      values=c('#40B0A6', 'coral'), 
                      limits = c(TRUE, FALSE),
                      labels=c('RP Acceptor', 'RP Non-Acceptor')),
  scale_fill_manual(name='Type', 
                    values=c('#b0e2de', '#ffccb9'), 
                    limits = c(TRUE, FALSE),
                    labels=c('RP Acceptor', 'RP Non-Acceptor')),
  theme(text = element_text(size = 12),
        legend.position = 'hide'),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE
  )
)

p1 + s
ggsave("plots/orc_ajt_acceptance_proficiency.png", width=6.5, height=2, dpi=600)
p2 + s
ggsave("plots/src_ajt_acceptance_proficiency.png", width=6.5, height=2, dpi=600)

md <- plot2 %>%
  select(study, group, participant, proficiency, acceptor) %>%
  group_by(study) %>%
  mutate(proficiency = scale(as.numeric(proficiency))) %>%
  ungroup() %>%
  filter(is.na(proficiency) == FALSE)

summary(md)

# apply deviation coding ...

contrasts(md$group) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

# view contrasts ...

contrasts(md$group)

# fit model ...

tic()
model <- glm(acceptor ~ group * proficiency, 
               data = filter(md, study == 'ORC'), 
               family = binomial)
summary(model)
toc()
beep()

# Parameter          Estimate Std. Error z value Pr(>|z|)   
# (Intercept)          0.1505     0.1734   0.868  0.38531   
# group2               0.1229     0.3468   0.354  0.72302   
# proficiency         -0.5285     0.1872  -2.823  0.00476 **
# group2:proficiency   0.3193     0.3745   0.853  0.39385  

# fit model ...

tic()
model <- glm(acceptor ~ group * proficiency, 
             data = filter(md, study == 'SRC'), 
             family = binomial)
summary(model)
toc()
beep()

# Parameter          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         0.86822    0.20445   4.247 2.17e-05 ***
# group2              1.34228    0.40890   3.283  0.00103 ** 
# proficiency         0.05472    0.20399   0.268  0.78852    
# group2:proficiency  0.25854    0.40798   0.634  0.52627

# ---------------------------------------------------------------------------- #
# plot: rt ####
# ---------------------------------------------------------------------------- #

# summarize to ajt data ...

temp1 <- ajt_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, 
                                '<br>(*n* = ', pairing_size, ')'))

# remove extreme values ...

temp2 <- temp1 %>%
  filter(rt <= 30000, rt >= 100)

# calculate percent excluded ...

check <- temp1 %>%
  mutate(keep = case_when(rt <= 30000 & rt >= 100 ~ TRUE,
                          TRUE ~ FALSE))

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, pairing) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# remove outliers more than 2 SDs from the mean ...

temp3 <- temp2 %>%
  group_by(study, pairing, dependency, environment) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  filter(rt <= upper & rt >= lower) %>%
  select(-upper, -lower)

# calculate percent excluded ...

check <- temp2 %>%
  group_by(study, pairing, dependency, environment, accuracy) %>%
  mutate(upper = mean(rt) + (2 * sd(rt)),
         lower = mean(rt) - (2 * sd(rt))) %>%
  ungroup() %>%
  mutate(keep = case_when(rt <= upper & rt >= lower ~ TRUE,
                          TRUE ~ FALSE)) %>%
  select(study, pairing, participant, dependency, environment, 
         accuracy, rt, upper, lower, keep)

check1 <- check %>%
  group_by(study) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

check2 <- check %>%
  group_by(study, pairing) %>%
  summarise(prop = length(which(keep == FALSE)) / n() * 100) %>%
  ungroup()

# define rt dataframe

ajt_rt <- temp3

# plot ...

plot <- ajt_rt %>%
  group_by(study, pairing, pairing_label, dependency, environment) %>%
  summarise(mean = mean(rt, na.rm=T),
            lcl = my.lcl(rt),
            ucl = my.ucl(rt)) %>%
  ungroup() %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = case_when(environment == 'Island' & pairing %in% c('KK', 'MM') ~ "‘Island’",
                                 TRUE ~ environment)) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island", "‘Island’"))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

design <- '
ABD
#CE
'

s <- list(
  aes(x=environment, y=mean/1000, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl/1000, ymax=ucl/1000), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Environment", 
                   #limits = c("short", "long", "island"), 
                   #labels = c("Short", "Long", "Island")
                   ),
  scale_y_continuous(name="RT (s)",
                     breaks = c(4, 6, 8, 10)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
  # theme(
  #   text = element_text(size = 12),
  #   plot.margin = margin(2, 2, 2, 2),
  #   legend.position = c(.15, .22),
  #   legend.margin = margin(-10, 0, 0, 0),
  #   axis.title.y = element_text(hjust = .88),
  #   strip.text = element_markdown()),
  # facet_manual(~pairing_label, 
  #              design = design, 
  #              axes = 'all', 
  #              remove_labels = FALSE) 
) 

p1 + s
ggsave("plots/orc_ajt_rt.png", width=6.5, height=2.25, dpi=600)

p2 + s
ggsave("plots/src_ajt_rt.png", width=6.5, height=2.25, dpi=600)

# table...

dt <- plot %>%
  select(study, pairing, dependency, environment, mean, lcl, ucl) %>%
  my_table() %>%
  csv_table()

dt %>% write_csv("tables/ajt_rt.csv")

# ---------------------------------------------------------------------------- #
# plot: logrt ####
# ---------------------------------------------------------------------------- #

# plot ...

plot <- ajt_rt %>%
  mutate(logrt = log(rt)) %>%
  group_by(study, pairing, pairing_label, dependency, environment) %>%
  summarise(mean = mean(logrt, na.rm=T),
            lcl = my.lcl(logrt),
            ucl = my.ucl(logrt)) %>%
  ungroup() %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = case_when(environment == 'Island' & pairing %in% c('KK', 'MM') ~ "‘Island’",
                                 TRUE ~ environment)) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island", "‘Island’"))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name="Environment", 
                   #limits = c("short", "long", "island"), 
                   #labels = c("Short", "Long", "Island")
  ),
  scale_y_continuous(name="Log RT", breaks = c(8, 8.5, 9)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
  # theme(
  #   text = element_text(size = 12),
  #   plot.margin = margin(2, 2, 2, 2),
  #   legend.position = c(.15, .22),
  #   legend.margin = margin(-10, 0, 0, 0),
  #   axis.title.y = element_text(hjust = .88),
  #   strip.text = element_markdown()),
  # facet_manual(~pairing_label, 
  #              design = design, 
  #              axes = 'all', 
  #              remove_labels = FALSE) 
) 

p1 + s + coord_cartesian(ylim = c(8, 9))
ggsave("plots/orc_ajt_logrt.png", width=6.5, height=2.25, dpi=600)

p2 + s + coord_cartesian(ylim = c(8, 9))
ggsave("plots/src_ajt_logrt.png", width=6.5, height=2.25, dpi=600)

# table...

dt <- plot %>%
  select(study, pairing, dependency, environment, mean, lcl, ucl) %>%
  my_table() %>%
  csv_table()

dt %>% write_csv("tables/ajt_logrt.csv")

# ---------------------------------------------------------------------------- #
# model: rt ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- ajt_rt %>%
  select(study, pairing, participant, item, dependency, environment, rt) %>%
  mutate_at(c('study', 'pairing', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$pairing) <- contr.treatment(5) - matrix(rep(1/5, 20), ncol = 4) 

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit ...

tic()
model2 <- lmer(rt ~ dependency * environment * pairing + 
                 (1 + dependency * environment | participant) +
                 (1 | item), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 62.37 sec

check_convergence(model2) # okay
check_singularity(model2) # warning (but random effects table looks okay) 

# set as final model ...

model <- model2

# save ...

model %>% write_rds('models/orc_ajt_rt_model.rds')
model <- read_rds('models/orc_ajt_rt_model.rds') 

# check performance ...

check_convergence(model) # convergence: okay
check_singularity(model, tolerance = 1e-9) # singularity: warning
check_collinearity(model) # independence: okay
check_heteroscedasticity(model) # homogeneity of variance: warning
check_outliers(model) # outliers: okay
check_autocorrelation(model) # autocorrelation: warning (normal)
MuMIn::r.squaredGLMM(model) # R2m = 0.1361187, R2c = 0.3384674
check_model(model) # linearity: okay; normality: a bit off
ggsave('model.pdf', width=10, height=15) 

# summary ...

summary(model)

# Parameter                         Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                        6511.17     105.00  222.56  62.009  < 2e-16 ***
# dependency2                         308.96      62.10  209.02   4.975 1.36e-06 ***
# environment2                        310.73      67.73  270.65   4.588 6.86e-06 ***
# environment3                        517.58      75.95  206.95   6.815 1.01e-10 ***
# pairing2                           2410.75     245.01  251.42   9.839  < 2e-16 ***
# pairing3                           2999.33     238.43  250.88  12.579  < 2e-16 ***
# pairing4                            323.45     244.68  250.07   1.322  0.18740    
# pairing5                           1467.81     238.37  250.59   6.158 2.92e-09 ***
# dependency2:environment2           -191.31     134.08  345.70  -1.427  0.15454    
# dependency2:environment3           -463.82     139.68  231.62  -3.321  0.00104 ** 
# dependency2:pairing2                107.40     180.99  530.69   0.593  0.55316    
# dependency2:pairing3                376.41     175.91  533.25   2.140  0.03282 *  
# dependency2:pairing4               1007.18     179.40  516.51   5.614 3.23e-08 ***
# dependency2:pairing5                511.02     175.31  522.99   2.915  0.00371 ** 
# environment2:pairing2               221.93     206.78  808.28   1.073  0.28346    
# environment3:pairing2                29.57     221.16  523.98   0.134  0.89370    
# environment2:pairing3               218.59     200.88  814.53   1.088  0.27685    
# environment3:pairing3               227.26     215.29  529.12   1.056  0.29163    
# environment2:pairing4               123.16     204.44  784.24   0.602  0.54708    
# environment3:pairing4               -77.11     219.62  512.68  -0.351  0.72564    
# environment2:pairing5               -57.30     200.20  797.58  -0.286  0.77480    
# environment3:pairing5              -137.03     214.69  519.83  -0.638  0.52360    
# dependency2:environment2:pairing2   -15.90     411.06 1034.02  -0.039  0.96916    
# dependency2:environment3:pairing2   109.85     420.37  662.13   0.261  0.79392    
# dependency2:environment2:pairing3   749.66     399.27 1041.38   1.878  0.06072 .  
# dependency2:environment3:pairing3  1010.48     409.27  670.59   2.469  0.01380 *  
# dependency2:environment2:pairing4  -744.46     406.87 1008.87  -1.830  0.06759 .  
# dependency2:environment3:pairing4  -643.78     417.30  648.25  -1.543  0.12339    
# dependency2:environment2:pairing5   -34.19     398.05 1021.68  -0.086  0.93157    
# dependency2:environment3:pairing5  -233.68     407.80  656.11  -0.573  0.56682   

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_ajt_rt_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rt_model.csv')
dt <- read_rds('tables/orc_ajt_rt_model.rds')
csv <- read_csv('tables/orc_ajt_rt_model.csv')

dt %>% kbl(format = 'markdown')

# |Parameter                                 |Estimate~($\hat{\beta}$) |$SE$   |$t$   |$p$    |
# |:-----------------------------------------|:------------------------|:------|:-----|:------|
# |(Intercept)                               |6511.17                  |105.00 |62.01 |< .001 |
# |Dependency[RP]                            |308.96                   |62.10  |4.98  |< .001 |
# |Environment[Long]                         |310.73                   |67.73  |4.59  |< .001 |
# |Environment[Island]                       |517.58                   |75.95  |6.81  |< .001 |
# |Pairing[KE]                               |2410.75                  |245.01 |9.84  |< .001 |
# |Pairing[ME]                               |2999.33                  |238.43 |12.58 |< .001 |
# |Pairing[KK]                               |323.45                   |244.68 |1.32  |0.187  |
# |Pairing[MM]                               |1467.81                  |238.37 |6.16  |< .001 |
# |Dependency[RP]:Environment[Long]          |−191.31                  |134.08 |−1.43 |0.155  |
# |Dependency[RP]:Environment[Island]        |−463.82                  |139.68 |−3.32 |0.001  |
# |Dependency[RP]:Pairing[KE]                |107.40                   |180.99 |0.59  |0.553  |
# |Dependency[RP]:Pairing[ME]                |376.41                   |175.91 |2.14  |0.033  |
# |Dependency[RP]:Pairing[KK]                |1007.18                  |179.40 |5.61  |< .001 |
# |Dependency[RP]:Pairing[MM]                |511.02                   |175.31 |2.91  |0.004  |
# |Environment[Long]:Pairing[KE]             |221.93                   |206.78 |1.07  |0.283  |
# |Environment[Island]:Pairing[KE]           |29.57                    |221.16 |0.13  |0.894  |
# |Environment[Long]:Pairing[ME]             |218.59                   |200.88 |1.09  |0.277  |
# |Environment[Island]:Pairing[ME]           |227.26                   |215.29 |1.06  |0.292  |
# |Environment[Long]:Pairing[KK]             |123.16                   |204.44 |0.60  |0.547  |
# |Environment[Island]:Pairing[KK]           |−77.11                   |219.62 |−0.35 |0.726  |
# |Environment[Long]:Pairing[MM]             |−57.30                   |200.20 |−0.29 |0.775  |
# |Environment[Island]:Pairing[MM]           |−137.03                  |214.69 |−0.64 |0.524  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |−15.90                   |411.06 |−0.04 |0.969  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |109.85                   |420.37 |0.26  |0.794  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |749.66                   |399.27 |1.88  |0.061  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |1010.48                  |409.27 |2.47  |0.014  |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |−744.46                  |406.87 |−1.83 |0.068  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−643.78                  |417.30 |−1.54 |0.123  |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |−34.19                   |398.05 |−0.09 |0.932  |
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |−233.68                  |407.80 |−0.57 |0.567  |

# joint tests ...

joint_tests(model)

# model term                     df1 df2 F.ratio p.value
# dependency                       1 Inf  24.753  <.0001 ***
# environment                      2 Inf  24.640  <.0001 ***
# pairing                          4 Inf 244.883  <.0001 ***
# dependency:environment           2 Inf   5.583  0.0038 ***
# dependency:pairing               4 Inf  10.190  <.0001 ***
# environment:pairing              8 Inf   0.670  0.7186
# dependency:environment:pairing   8 Inf   2.516  0.0099 **

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate  SE  df z.ratio p.value
# pronoun - gap short       EE        140.13 197 Inf   0.712  0.9999
# pronoun - gap long        EE        -42.20 195 Inf  -0.217  1.0000
# pronoun - gap island      EE       -372.27 196 Inf  -1.903  0.5837
# pronoun - gap short       KE        216.21 228 Inf   0.950  0.9981
# pronoun - gap long        KE         17.98 228 Inf   0.079  1.0000
# pronoun - gap island      KE       -186.33 226 Inf  -0.823  0.9996
# pronoun - gap short       ME        -70.17 216 Inf  -0.324  1.0000
# pronoun - gap long        ME        497.16 216 Inf   2.305  0.2740
# pronoun - gap island      ME        427.91 216 Inf   1.978  0.5195
# pronoun - gap short       KK       1610.05 225 Inf   7.161  <.0001 ***
# pronoun - gap long        KK        683.26 222 Inf   3.075  0.0311 *
# pronoun - gap island      KK        453.88 223 Inf   2.031  0.4752
# pronoun - gap short       MM        740.44 216 Inf   3.424  0.0092 **
# pronoun - gap long        MM        523.92 214 Inf   2.453  0.1922
# pronoun - gap island      MM         -5.63 214 Inf  -0.026  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_rt_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rt_posthoc.csv') 
dt <- read_rds('tables/orc_ajt_rt_posthoc.rds')
csv <- read_csv('tables/orc_ajt_rt_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit ...

tic()
model3 <- lmer(rt ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant) +
                 (1 + environment | item), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

# 35.18 sec

check_convergence(model3) # okay
check_singularity(model3) # warning (but random effects table looks okay)

# set as final model ...

model <- model3

# save ...

model %>% write_rds('models/src_ajt_rt_model.rds')
model <- read_rds('models/src_ajt_rt_model.rds')

# check performance ...

check_convergence(model) # convergence: okay
check_singularity(model, tolerance = 1e-9) # singularity: warning
check_collinearity(model) # independence: okay
check_heteroscedasticity(model) # homogeneity of variance: warning
check_outliers(model) # outliers: okay
check_autocorrelation(model) # autocorrelation: warning (normal)
MuMIn::r.squaredGLMM(model) # R2m = 0.2035956 , R2c = 0.4005973
check_model(model) # linearity: okay; normality: a bit off
ggsave('model.pdf', width=10, height=15) 

# summary ...

summary(model)
  
# Parameter                         Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                        6257.97     112.79   199.47  55.482  < 2e-16 ***
# dependency2                         658.54      64.42   227.64  10.223  < 2e-16 ***
# environment2                        486.02      74.18    76.66   6.552 5.88e-09 ***
# environment3                        869.53      84.80   141.76  10.254  < 2e-16 ***
# pairing2                           2395.66     274.19   214.50   8.737 6.91e-16 ***
# pairing3                           3707.39     268.10   215.23  13.829  < 2e-16 ***
# pairing4                            222.02     273.95   213.78   0.810 0.418595    
# pairing5                           1947.86     267.71   214.02   7.276 6.42e-12 ***
# dependency2:environment2          -1206.35     133.28  8820.26  -9.051  < 2e-16 ***
# dependency2:environment3          -1348.70     133.33  8820.39 -10.116  < 2e-16 ***
# dependency2:pairing2                185.72     197.15   535.75   0.942 0.346618    
# dependency2:pairing3                368.73     193.42   542.34   1.906 0.057135 .  
# dependency2:pairing4                894.05     196.08   527.35   4.560 6.37e-06 ***
# dependency2:pairing5                836.52     191.88   529.87   4.360 1.57e-05 ***
# environment2:pairing2               330.47     228.27   841.94   1.448 0.148065    
# environment3:pairing2               385.34     249.74   431.72   1.543 0.123570    
# environment2:pairing3               -19.07     224.21   856.39  -0.085 0.932246    
# environment3:pairing3               571.85     245.27   437.77   2.331 0.020180 *  
# environment2:pairing4               160.62     226.64   826.69   0.709 0.478704    
# environment3:pairing4               185.52     248.80   427.50   0.746 0.456272    
# environment2:pairing5               -23.48     221.86   831.77  -0.106 0.915734    
# environment3:pairing5               136.22     243.29   428.48   0.560 0.575835    
# dependency2:environment2:pairing2  -144.98     435.60  8849.28  -0.333 0.739264    
# dependency2:environment3:pairing2 -1633.88     434.94  8834.82  -3.757 0.000173 ***
# dependency2:environment2:pairing3  -647.40     427.54  8823.56  -1.514 0.130004    
# dependency2:environment3:pairing3 -1571.79     427.05  8820.48  -3.681 0.000234 ***
# dependency2:environment2:pairing4  -303.16     432.41  8839.97  -0.701 0.483266    
# dependency2:environment3:pairing4 -1177.99     432.78  8823.11  -2.722 0.006503 ** 
# dependency2:environment2:pairing5 -1967.62     422.95  8823.82  -4.652 3.33e-06 ***
# dependency2:environment3:pairing5 -3387.46     423.28  8816.36  -8.003 1.37e-15 ***

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_ajt_rt_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rt_model.csv')
dt <- read_rds('tables/src_ajt_rt_model.rds')
csv <- read_csv('tables/src_ajt_rt_model.csv')

dt %>% kbl(format = 'markdown')

# |Parameter                                 |Estimate~($\hat{\beta}$) |$SE$   |$t$    |$p$    |
# |:-----------------------------------------|:------------------------|:------|:------|:------|
# |(Intercept)                               |6257.97                  |112.79 |55.48  |< .001 |
# |Dependency[RP]                            |658.54                   |64.42  |10.22  |< .001 |
# |Environment[Long]                         |486.02                   |74.18  |6.55   |< .001 |
# |Environment[Island]                       |869.53                   |84.80  |10.25  |< .001 |
# |Pairing[KE]                               |2395.66                  |274.19 |8.74   |< .001 |
# |Pairing[ME]                               |3707.39                  |268.10 |13.83  |< .001 |
# |Pairing[KK]                               |222.02                   |273.95 |0.81   |0.419  |
# |Pairing[MM]                               |1947.86                  |267.71 |7.28   |< .001 |
# |Dependency[RP]:Environment[Long]          |−1206.35                 |133.28 |−9.05  |< .001 |
# |Dependency[RP]:Environment[Island]        |−1348.70                 |133.33 |−10.12 |< .001 |
# |Dependency[RP]:Pairing[KE]                |185.72                   |197.15 |0.94   |0.347  |
# |Dependency[RP]:Pairing[ME]                |368.73                   |193.42 |1.91   |0.057  |
# |Dependency[RP]:Pairing[KK]                |894.05                   |196.08 |4.56   |< .001 |
# |Dependency[RP]:Pairing[MM]                |836.52                   |191.88 |4.36   |< .001 |
# |Environment[Long]:Pairing[KE]             |330.47                   |228.27 |1.45   |0.148  |
# |Environment[Island]:Pairing[KE]           |385.34                   |249.74 |1.54   |0.124  |
# |Environment[Long]:Pairing[ME]             |−19.07                   |224.21 |−0.09  |0.932  |
# |Environment[Island]:Pairing[ME]           |571.85                   |245.27 |2.33   |0.020  |
# |Environment[Long]:Pairing[KK]             |160.62                   |226.64 |0.71   |0.479  |
# |Environment[Island]:Pairing[KK]           |185.52                   |248.80 |0.75   |0.456  |
# |Environment[Long]:Pairing[MM]             |−23.48                   |221.86 |−0.11  |0.916  |
# |Environment[Island]:Pairing[MM]           |136.22                   |243.29 |0.56   |0.576  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |−144.98                  |435.60 |−0.33  |0.739  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |−1633.88                 |434.94 |−3.76  |< .001 |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |−647.40                  |427.54 |−1.51  |0.130  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |−1571.79                 |427.05 |−3.68  |< .001 |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |−303.16                  |432.41 |−0.70  |0.483  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−1177.99                 |432.78 |−2.72  |0.007  |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |−1967.62                 |422.95 |−4.65  |< .001 |
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |−3387.46                 |423.28 |−8.00  |< .001 |

# joint tests ...

joint_tests(model)

# model term                     df1 df2 F.ratio p.value
# dependency                       1 Inf 104.513  <.0001 ***
# environment                      2 Inf  53.259  <.0001 ***
# pairing                          4 Inf 295.489  <.0001 ***
# dependency:environment           2 Inf  61.992  <.0001 ***
# dependency:pairing               4 Inf   9.445  <.0001 ***
# environment:pairing              8 Inf   1.818  0.0686 .
# dependency:environment:pairing   8 Inf   9.377  <.0001 ***

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate  SE  df z.ratio p.value
# pronoun - gap short       EE         330.9 228 Inf   1.449  0.9070
# pronoun - gap long        EE        -262.8 231 Inf  -1.137  0.9877
# pronoun - gap island      EE         536.5 231 Inf   2.327  0.2602
# pronoun - gap short       KE        1109.6 221 Inf   5.019  <.0001 ***
# pronoun - gap long        KE         370.9 223 Inf   1.662  0.7795
# pronoun - gap island      KE        -318.8 222 Inf  -1.435  0.9132
# pronoun - gap short       ME        1439.4 212 Inf   6.785  <.0001 ***
# pronoun - gap long        ME         198.3 215 Inf   0.924  0.9986
# pronoun - gap island      ME          73.1 214 Inf   0.342  1.0000
# pronoun - gap short       KK        1718.7 218 Inf   7.867  <.0001 ***
# pronoun - gap long        KK         821.8 219 Inf   3.756  0.0026 ***
# pronoun - gap island      KK         746.2 220 Inf   3.385  0.0106 *
# pronoun - gap short       MM        2952.5 208 Inf  14.220  <.0001 ***
# pronoun - gap long        MM         391.1 210 Inf   1.866  0.6156
# pronoun - gap island      MM        -229.5 211 Inf  -1.090  0.9919

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_rt_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rt_posthoc.csv') 
dt <- read_rds('tables/src_ajt_rt_posthoc.rds')
csv <- read_csv('tables/src_ajt_rt_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: log rt ####
# ---------------------------------------------------------------------------- #

# use md dataframe from raw rt analysis

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit ... 

tic()
model1 <- lmer(log(rt) ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant) +
                 (1 + dependency + environment | item), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 72 sec

check_convergence(model1) # okay
check_singularity(model1) # warning (but random effects table looks okay)

# set as final model ...

model <- model1

# save ...

# model %>% write_rds('models/orc_ajt_logrt_model.rds')
model <- read_rds('models/orc_ajt_logrt_model.rds')

# check performance ...

check_convergence(model) # convergence: okay
check_singularity(model, tolerance = 1e-9) # singularity: warning
check_collinearity(model) # independence: okay
check_heteroscedasticity(model) # homogeneity of variance: okay
check_outliers(model) # outliers: okay
check_autocorrelation(model) # autocorrelation: warning (normal)
MuMIn::r.squaredGLMM(model) # R2m = 0.1287457, R2c = 0.3750859
check_model(model) # linearity: okay; normality: a bit off
ggsave('model.pdf', width=10, height=15) 

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_ajt_logrt_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_logrt_model.csv')
dt <- read_rds('tables/orc_ajt_logrt_model.rds')
csv <- read_csv('tables/orc_ajt_logrt_model.csv')

dt %>% kbl(format = 'markdown')

# |Parameter                                 |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:-----------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                               |8.66                     |0.02 |500.93 |< .001 |
# |Dependency[RP]                            |0.05                     |0.01 |5.37   |< .001 |
# |Environment[Long]                         |0.04                     |0.01 |4.17   |< .001 |
# |Environment[Island]                       |0.07                     |0.01 |5.91   |< .001 |
# |Pairing[KE]                               |0.34                     |0.04 |8.50   |< .001 |
# |Pairing[ME]                               |0.43                     |0.04 |11.02  |< .001 |
# |Pairing[KK]                               |0.03                     |0.04 |0.81   |0.417  |
# |Pairing[MM]                               |0.22                     |0.04 |5.65   |< .001 |
# |Dependency[RP]:Environment[Long]          |−0.03                    |0.02 |−1.86  |0.063  |
# |Dependency[RP]:Environment[Island]        |−0.07                    |0.02 |−3.99  |< .001 |
# |Dependency[RP]:Pairing[KE]                |−0.01                    |0.03 |−0.36  |0.722  |
# |Dependency[RP]:Pairing[ME]                |0.03                     |0.03 |1.15   |0.250  |
# |Dependency[RP]:Pairing[KK]                |0.14                     |0.03 |5.29   |< .001 |
# |Dependency[RP]:Pairing[MM]                |0.05                     |0.03 |1.77   |0.077  |
# |Environment[Long]:Pairing[KE]             |0.02                     |0.03 |0.62   |0.536  |
# |Environment[Island]:Pairing[KE]           |−0.01                    |0.03 |−0.42  |0.675  |
# |Environment[Long]:Pairing[ME]             |0.01                     |0.03 |0.25   |0.802  |
# |Environment[Island]:Pairing[ME]           |−0.01                    |0.03 |−0.25  |0.804  |
# |Environment[Long]:Pairing[KK]             |0.02                     |0.03 |0.65   |0.516  |
# |Environment[Island]:Pairing[KK]           |−0.01                    |0.03 |−0.36  |0.720  |
# |Environment[Long]:Pairing[MM]             |−0.02                    |0.03 |−0.68  |0.498  |
# |Environment[Island]:Pairing[MM]           |−0.03                    |0.03 |−0.98  |0.328  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |0.02                     |0.06 |0.35   |0.728  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |0.07                     |0.06 |1.26   |0.208  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |0.10                     |0.06 |1.78   |0.074  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |0.15                     |0.06 |2.66   |0.008  |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |−0.13                    |0.06 |−2.24  |0.025  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−0.11                    |0.06 |−1.96  |0.050  |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |−0.01                    |0.06 |−0.11  |0.913  |
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |−0.01                    |0.06 |−0.16  |0.869  |

# joint tests ...

joint_tests(model)

# model term                     df1 df2 F.ratio p.value
# dependency                       1 Inf  28.865  <.0001 ***
# environment                      2 Inf  18.516  <.0001 ***
# pairing                          4 Inf 234.250  <.0001 ***
# dependency:environment           2 Inf   7.957  0.0004 ***
# dependency:pairing               4 Inf  10.939  <.0001 ***
# environment:pairing              8 Inf   0.402  0.9200
# dependency:environment:pairing   8 Inf   3.163  0.0014 **

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing  estimate     SE  df z.ratio p.value
# pronoun - gap short       EE       0.053374 0.0282 Inf   1.892  0.5912
# pronoun - gap long        EE       0.021228 0.0283 Inf   0.749  0.9999
# pronoun - gap island      EE      -0.041943 0.0284 Inf  -1.475  0.8938
# pronoun - gap short       KE       0.012488 0.0326 Inf   0.383  1.0000
# pronoun - gap long        KE       0.000644 0.0331 Inf   0.019  1.0000
# pronoun - gap island      KE      -0.009298 0.0329 Inf  -0.283  1.0000
# pronoun - gap short       ME      -0.000633 0.0310 Inf  -0.020  1.0000
# pronoun - gap long        ME       0.068483 0.0314 Inf   2.184  0.3545
# pronoun - gap island      ME       0.055294 0.0314 Inf   1.759  0.7030
# pronoun - gap short       KK       0.275822 0.0322 Inf   8.569  <.0001 ***
# pronoun - gap long        KK       0.114071 0.0323 Inf   3.530  0.0062 **
# pronoun - gap island      KK       0.067011 0.0325 Inf   2.064  0.4467
# pronoun - gap short       MM       0.104796 0.0310 Inf   3.382  0.0107 *
# pronoun - gap long        MM       0.066477 0.0311 Inf   2.140  0.3866
# pronoun - gap island      MM       0.000149 0.0311 Inf   0.005  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_logrt_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_logrt_posthoc.csv') 
dt <- read_rds('tables/orc_ajt_logrt_posthoc.rds')
csv <- read_csv('tables/orc_ajt_logrt_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit ... 

tic()
model5 <- lmer(log(rt) ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant) +
                 (1 + pairing | item), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model5)
toc()
beep()

# 802.1 sec elapsed
# singular fit

check_convergence(model5) # warning, but lme4 doesn't throw a warning
check_singularity(model5) # warning, but random-effects table looks okay

# save ...

# model5 %>% write_rds('models/src_ajt_logrt_model.rds')
model <- read_rds('models/src_ajt_logrt_model.rds')

# check performance ...

check_singularity(model, tolerance = 1e-9) # singularity: okay
check_collinearity(model) # independence: okay
check_heteroscedasticity(model) # homogeneity of variance: okay
check_outliers(model) # outliers: okay
check_autocorrelation(model) # autocorrelation: warning (normal)
MuMIn::r.squaredGLMM(model) # R2m = 0.204686, R2c = 0.4481549
check_model(model) # linearity: okay; normality: okay
ggsave('model.pdf', width=10, height=15) 

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_ajt_logrt_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_logrt_model.csv')
dt <- read_rds('tables/src_ajt_logrt_model.rds')
csv <- read_csv('tables/src_ajt_logrt_model.csv')

dt %>% kbl(format = 'markdown')

# |Parameter                                 |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:-----------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                               |8.60                     |0.02 |449.90 |< .001 |
# |Dependency[RP]                            |0.10                     |0.01 |10.22  |< .001 |
# |Environment[Long]                         |0.08                     |0.01 |6.93   |< .001 |
# |Environment[Island]                       |0.14                     |0.01 |10.84  |< .001 |
# |Pairing[KE]                               |0.37                     |0.05 |7.82   |< .001 |
# |Pairing[ME]                               |0.55                     |0.05 |11.79  |< .001 |
# |Pairing[KK]                               |0.00                     |0.05 |−0.10  |0.923  |
# |Pairing[MM]                               |0.31                     |0.05 |6.48   |< .001 |
# |Dependency[RP]:Environment[Long]          |−0.22                    |0.02 |−10.84 |< .001 |
# |Dependency[RP]:Environment[Island]        |−0.22                    |0.02 |−11.24 |< .001 |
# |Dependency[RP]:Pairing[KE]                |0.01                     |0.03 |0.49   |0.623  |
# |Dependency[RP]:Pairing[ME]                |0.02                     |0.03 |0.73   |0.467  |
# |Dependency[RP]:Pairing[KK]                |0.17                     |0.03 |5.74   |< .001 |
# |Dependency[RP]:Pairing[MM]                |0.12                     |0.03 |4.07   |< .001 |
# |Environment[Long]:Pairing[KE]             |0.02                     |0.03 |0.58   |0.564  |
# |Environment[Island]:Pairing[KE]           |0.02                     |0.04 |0.53   |0.593  |
# |Environment[Long]:Pairing[ME]             |−0.03                    |0.03 |−0.80  |0.424  |
# |Environment[Island]:Pairing[ME]           |0.02                     |0.04 |0.44   |0.662  |
# |Environment[Long]:Pairing[KK]             |0.03                     |0.03 |0.77   |0.444  |
# |Environment[Island]:Pairing[KK]           |0.04                     |0.04 |1.03   |0.302  |
# |Environment[Long]:Pairing[MM]             |−0.01                    |0.03 |−0.32  |0.748  |
# |Environment[Island]:Pairing[MM]           |0.02                     |0.04 |0.57   |0.568  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |−0.01                    |0.07 |−0.16  |0.875  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |−0.23                    |0.07 |−3.59  |< .001 |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |−0.01                    |0.06 |−0.18  |0.856  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |−0.20                    |0.06 |−3.18  |0.001  |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |−0.09                    |0.06 |−1.35  |0.176  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−0.27                    |0.06 |−4.12  |< .001 |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |−0.26                    |0.06 |−4.11  |< .001 |
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |−0.48                    |0.06 |−7.54  |< .001 |

# joint tests ...

joint_tests(model)

# model term                     df1 df2 F.ratio p.value
# dependency                       1 Inf 104.446  <.0001
# environment                      2 Inf  58.877  <.0001
# pairing                          4 Inf 246.344  <.0001
# dependency:environment           2 Inf  81.585  <.0001
# dependency:pairing               4 Inf  16.149  <.0001
# environment:pairing              8 Inf   0.632  0.7512
# dependency:environment:pairing   8 Inf   8.342  <.0001

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate     SE  df z.ratio p.value
# pronoun - gap short       EE        0.0801 0.0345 Inf   2.320  0.2638
# pronoun - gap long        EE       -0.0623 0.0349 Inf  -1.785  0.6820
# pronoun - gap island      EE        0.0920 0.0348 Inf   2.641  0.1166
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# pronoun - gap short       KE        0.1762 0.0334 Inf   5.281  <.0001 *
# pronoun - gap long        KE        0.0236 0.0337 Inf   0.699  0.9999
# pronoun - gap island      KE       -0.0454 0.0336 Inf  -1.354  0.9432
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# pronoun - gap short       ME        0.1734 0.0321 Inf   5.408  <.0001 *
# pronoun - gap long        ME        0.0194 0.0324 Inf   0.597  1.0000
# pronoun - gap island      ME       -0.0182 0.0323 Inf  -0.563  1.0000
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# pronoun - gap short       KK        0.3706 0.0330 Inf  11.232  <.0001 *
# pronoun - gap long        KK        0.1406 0.0331 Inf   4.253  0.0003 *
# pronoun - gap island      KK        0.1160 0.0333 Inf   3.482  0.0074 *
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# pronoun - gap short       MM        0.4458 0.0314 Inf  14.203  <.0001 *
# pronoun - gap long        MM        0.0432 0.0317 Inf   1.362  0.9406
# pronoun - gap island      MM       -0.0203 0.0318 Inf  -0.637  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_logrt_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_logrt_posthoc.csv') 
dt <- read_rds('tables/src_ajt_logrt_posthoc.rds')
csv <- read_csv('tables/src_ajt_logrt_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: rt proficiency ####
# ---------------------------------------------------------------------------- #

# prep dataframe

md <- ajt_rt %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(!condition %in% c('grammatical', 'ungrammatical')) %>%
  select(study, pairing, participant, item, dependency, environment, response, rt, proficiency) %>%
  mutate_at(c('study', 'pairing', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(response = ordered(response))

# remove participants with no proficiency scores ...

md <- md %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(proficiency = as.numeric(proficiency))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol = 2)
contrasts(md$pairing) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# plot ...

plot <- md %>%
  mutate(group = case_when(pairing == 'KE' ~ 'KLE', 
                           pairing == 'ME' ~ 'MLE')) %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in English<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, participant, proficiency, dependency, environment) %>%
  summarise(rt = mean(rt, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = rt/1000, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='RT (s)'),
  #coord_cartesian(ylim = c(0.9, 6.1)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(4, 10))
ggsave('plots/orc_ajt_rt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(4, 10))
ggsave('plots/src_ajt_rt_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: logrt proficiency ####
# ---------------------------------------------------------------------------- #

s <- list(
  aes(x = proficiency, y = log(rt), col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='Log RT', breaks = c(8.5, 9, 9.5)),
  #coord_cartesian(ylim = c(0.9, 6.1)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(8.4, 9.6))
ggsave('plots/orc_ajt_logrt_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(8.4, 9.6))
ggsave('plots/src_ajt_logrt_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# model: rt proficiency ####
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic() 
model3 <- lmer(rt ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 + environment | item),
               data = filter(md, study == 'ORC'), 
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

# 16 sec

check_convergence(model3) # okay
check_singularity(model3) # warning (but random effects table looks okay)

# set as final model ...

model <- model3

# save ...

model %>% write_rds('models/orc_ajt_rt_proficiency_model.rds')
model <- read_rds('models/orc_ajt_rt_proficiency_model.rds')

# model summary ...

summary(model)

# Parameter                                     Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                   7821.412    179.864  137.495  43.485  < 2e-16 ***
# dependency2                                    133.899    118.195  140.163   1.133  0.25921    
# environment2                                   448.045    133.003   46.511   3.369  0.00153 ** 
# environment3                                   656.792    135.822   56.644   4.836 1.06e-05 ***
# pairing2                                       516.110    357.562  137.987   1.443  0.15117    
# proficiency                                    205.831    180.223  137.134   1.142  0.25541    
# dependency2:environment2                       189.916    252.015 3467.445   0.754  0.45115    
# dependency2:environment3                        42.321    251.858 3469.668   0.168  0.86657    
# dependency2:pairing2                           276.571    236.359  139.997   1.170  0.24394    
# environment2:pairing2                           -6.522    258.965  182.367  -0.025  0.97994    
# environment3:pairing2                          189.614    257.311  864.926   0.737  0.46138    
# dependency2:proficiency                         98.194    118.250  136.277   0.830  0.40777    
# environment2:proficiency                       129.250    129.268  176.838   1.000  0.31875    
# environment3:proficiency                       255.938    127.783  840.731   2.003  0.04551 *  
# pairing2:proficiency                           464.769    360.445  137.135   1.289  0.19942    
# dependency2:environment2:pairing2              754.387    503.915 3474.964   1.497  0.13447    
# dependency2:environment3:pairing2              927.853    503.535 3474.948   1.843  0.06546 .  
# dependency2:environment2:proficiency            24.365    251.799 3381.716   0.097  0.92292    
# dependency2:environment3:proficiency          -294.423    250.650 3432.355  -1.175  0.24022    
# dependency2:pairing2:proficiency               -64.136    236.272  135.886  -0.271  0.78646    
# environment2:pairing2:proficiency              122.055    258.564  177.000   0.472  0.63747    
# environment3:pairing2:proficiency             -253.466    255.607  841.602  -0.992  0.32167    
# dependency2:environment2:pairing2:proficiency  180.345    503.367 3463.588   0.358  0.72016    
# dependency2:environment3:pairing2:proficiency -450.854    500.743 3465.993  -0.900  0.36799 

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_rt_proficiency_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rt_proficiency_model.csv') 
dt <- read_rds('tables/orc_ajt_rt_proficiency_model.rds')

# joint tests ...

joint_tests(model)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 Inf   1.283  0.2574
# environment                                  2 Inf  12.426  <.0001 ***
# pairing                                      1 Inf   2.082  0.1490
# proficiency                                  1 Inf   1.304  0.2534
# dependency:environment                       2 Inf   0.311  0.7327
# dependency:pairing                           1 Inf   1.369  0.2419
# dependency:proficiency                       1 Inf   0.690  0.4063
# environment:pairing                          2 Inf   0.356  0.7004
# environment:proficiency                      2 Inf   2.010  0.1340
# pairing:proficiency                          1 Inf   1.663  0.1972
# dependency:environment:pairing               2 Inf   1.927  0.1456
# dependency:environment:proficiency           2 Inf   0.991  0.3711
# dependency:pairing:proficiency               1 Inf   0.074  0.7860
# environment:pairing:proficiency              2 Inf   1.039  0.3537
# dependency:environment:pairing:proficiency   2 Inf   0.831  0.4358

pairwise <- model %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment pairing estimate  SE  df z.ratio p.value
# pronoun - gap short       KE           175 277 Inf   0.632  0.9887
# pronoun - gap long        KE           109 283 Inf   0.387  0.9992
# pronoun - gap island      KE           106 279 Inf   0.380  0.9993
# pronoun - gap short       ME           201 247 Inf   0.816  0.9589
# pronoun - gap long        ME           316 251 Inf   1.258  0.7520
# pronoun - gap island      ME          -319 250 Inf  -1.275  0.7398

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_rt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_rt_proficiency_posthoc.csv') 
dt <- read_rds('tables/orc_ajt_rt_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

tic() 
model2 <- lmer(rt ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 | item),
               data = filter(md, study == 'SRC'), 
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 4 sec

check_convergence(model2) # okay
check_singularity(model2) # warning (but random effects table looks okay)

# set as final model ...

model <- model2

# save ...

model %>% write_rds('models/src_ajt_rt_proficiency_model.rds')
model <- read_rds('models/src_ajt_rt_proficiency_model.rds')

# model summary ...

summary(model)

# Parameter                                     Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                    7677.92     199.59   133.94  38.468  < 2e-16 ***
# dependency2                                     502.41     112.49   144.67   4.466 1.59e-05 ***
# environment2                                    530.80     127.52   251.35   4.163 4.33e-05 ***
# environment3                                   1053.85     158.29   130.70   6.658 6.98e-10 ***
# pairing2                                       1266.57     390.95   129.93   3.240  0.00152 ** 
# proficiency                                     204.15     196.63   129.97   1.038  0.30109    
# dependency2:environment2                      -1001.46     245.88  3366.39  -4.073 4.75e-05 ***
# dependency2:environment3                      -1395.33     244.79  3355.79  -5.700 1.30e-08 ***
# dependency2:pairing2                            270.01     224.98   144.67   1.200  0.23204    
# environment2:pairing2                          -368.44     255.17   252.00  -1.444  0.15001    
# environment3:pairing2                           203.33     316.63   130.79   0.642  0.52187    
# dependency2:proficiency                         177.34     112.81   142.39   1.572  0.11816    
# environment2:proficiency                         47.42     127.22   243.30   0.373  0.70966    
# environment3:proficiency                        197.13     159.50   130.71   1.236  0.21871    
# pairing2:proficiency                             67.30     393.26   129.97   0.171  0.86439    
# dependency2:environment2:pairing2              -349.10     491.13  3363.15  -0.711  0.47725    
# dependency2:environment3:pairing2               105.26     489.73  3361.34   0.215  0.82984    
# dependency2:environment2:proficiency            275.61     244.41  3333.55   1.128  0.25954    
# dependency2:environment3:proficiency             26.09     246.31  3358.27   0.106  0.91564    
# dependency2:pairing2:proficiency                149.58     225.65   142.42   0.663  0.50846    
# environment2:pairing2:proficiency              -155.82     254.43   243.23  -0.612  0.54083    
# environment3:pairing2:proficiency               -11.12     319.00   130.73  -0.035  0.97223    
# dependency2:environment2:pairing2:proficiency  -242.57     489.02  3340.55  -0.496  0.61990    
# dependency2:environment3:pairing2:proficiency  -257.96     492.60  3358.44  -0.524  0.60054

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_rt_proficiency_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rt_proficiency_model.csv') 
dt <- read_rds('tables/src_ajt_rt_proficiency_model.rds')

# joint tests ...

joint_tests(model)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 Inf  20.055  <.0001 ***
# environment                                  2 Inf  22.712  <.0001 ***
# pairing                                      1 Inf  10.507  0.0012 **
# proficiency                                  1 Inf   1.078  0.2992
# dependency:environment                       2 Inf  17.271  <.0001 ***
# dependency:pairing                           1 Inf   1.452  0.2282
# dependency:proficiency                       1 Inf   2.471  0.1159
# environment:pairing                          2 Inf   2.280  0.1023
# environment:proficiency                      2 Inf   0.802  0.4484
# pairing:proficiency                          1 Inf   0.029  0.8641
# dependency:environment:pairing               2 Inf   0.468  0.6265
# dependency:environment:proficiency           2 Inf   0.771  0.4626
# dependency:pairing:proficiency               1 Inf   0.439  0.5074
# environment:pairing:proficiency              2 Inf   0.235  0.7902
# dependency:environment:pairing:proficiency   2 Inf   0.174  0.8400

pairwise <- model %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment pairing estimate  SE  df z.ratio p.value
# pronoun - gap short       KE         -81.4 263 Inf  -0.309  0.9998
# pronoun - gap long        KE         315.4 264 Inf   1.193  0.7948
# pronoun - gap island      KE          73.6 268 Inf   0.274  0.9999
# pronoun - gap short       ME         235.0 246 Inf   0.956  0.9157
# pronoun - gap long        ME         389.3 246 Inf   1.581  0.5141
# pronoun - gap island      ME         132.1 249 Inf   0.530  0.9956

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_rt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_rt_proficiency_posthoc.csv') 
dt <- read_rds('tables/src_ajt_rt_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: logrt proficiency ####
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic() 
model1 <- lmer(log(rt) ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 | item),
               data = filter(md, study == 'ORC'), 
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 4 sec
# no warnings

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/orc_ajt_logrt_proficiency_model.rds')
model <- read_rds('models/orc_ajt_logrt_proficiency_model.rds')

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_logrt_proficiency_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_logrt_proficiency_model.csv') 
dt <- read_rds('tables/orc_ajt_logrt_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                             |Est.  |$SE$ |df      |$t$    |$p$    |
# |:-----------------------------------------------------|:-----|:----|:-------|:------|:------|
# |(Intercept)                                           |8.84  |0.03 |136.07  |353.46 |< .001 | ***
# |Dependency[RP]                                        |0.02  |0.02 |142.40  |1.16   |0.247  |
# |Environment[Long]                                     |0.05  |0.02 |134.71  |3.27   |0.001  | **
# |Environment[Island]                                   |0.07  |0.02 |137.20  |4.51   |< .001 | ***
# |Pairing[ME]                                           |0.08  |0.05 |139.26  |1.56   |0.121  |
# |Proficiency                                           |0.03  |0.03 |138.60  |1.36   |0.178  |
# |Dependency[RP]:Environment[Long]                      |0.03  |0.03 |3380.46 |0.79   |0.431  |
# |Dependency[RP]:Environment[Island]                    |0.01  |0.03 |3379.87 |0.44   |0.659  |
# |Dependency[RP]:Pairing[ME]                            |0.04  |0.03 |142.36  |1.31   |0.191  |
# |Environment[Long]:Pairing[ME]                         |−0.01 |0.03 |134.30  |−0.29  |0.772  |
# |Environment[Island]:Pairing[ME]                       |0.00  |0.03 |136.62  |0.14   |0.891  |
# |Dependency[RP]:Proficiency                            |0.01  |0.02 |138.58  |0.52   |0.607  |
# |Environment[Long]:Proficiency                         |0.01  |0.02 |130.87  |0.71   |0.480  |
# |Environment[Island]:Proficiency                       |0.04  |0.02 |131.21  |2.24   |0.027  | *
# |Pairing[ME]:Proficiency                               |0.07  |0.05 |138.60  |1.47   |0.145  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]               |0.08  |0.06 |3385.24 |1.23   |0.217  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME]             |0.08  |0.06 |3377.49 |1.25   |0.212  |
# |Dependency[RP]:Envir.[Long]:Proficiency               |0.00  |0.03 |3375.87 |0.00   |> .999 |
# |Dependency[RP]:Envir.[Island]:Proficiency             |−0.04 |0.03 |3333.58 |−1.31  |0.191  |
# |Dependency[RP]:Pairing[ME]:Proficiency                |−0.01 |0.03 |138.40  |−0.19  |0.846  |
# |Environment[Long]:Pairing[ME]:Proficiency             |0.00  |0.03 |130.84  |−0.03  |0.972  |
# |Environment[Island]:Pairing[ME]:Proficiency           |−0.04 |0.03 |131.18  |−1.35  |0.180  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]:Proficiency   |0.01  |0.06 |3369.18 |0.10   |0.919  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME]:Proficiency |−0.04 |0.06 |3370.20 |−0.70  |0.481  |

# joint tests ...

joint_tests(model)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 Inf   1.351  0.2451
# environment                                  2 Inf  10.802  <.0001 ***
# pairing                                      1 Inf   2.435  0.1187
# proficiency                                  1 Inf   1.836  0.1754
# dependency:environment                       2 Inf   0.312  0.7320
# dependency:pairing                           1 Inf   1.724  0.1891
# dependency:proficiency                       1 Inf   0.266  0.6057
# environment:pairing                          2 Inf   0.097  0.9072
# environment:proficiency                      2 Inf   2.638  0.0715
# pairing:proficiency                          1 Inf   2.147  0.1429
# dependency:environment:pairing               2 Inf   1.031  0.3567
# dependency:environment:proficiency           2 Inf   1.131  0.3226
# dependency:pairing:proficiency               1 Inf   0.038  0.8457
# environment:pairing:proficiency              2 Inf   1.197  0.3022
# dependency:environment:pairing:proficiency   2 Inf   0.381  0.6829

# post-hoc tests ...

pairwise <- model %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing  estimate     SE  df z.ratio p.value
# pronoun - gap short       KE       0.018679 0.0359 Inf   0.520  0.9959
# pronoun - gap long        KE       0.015472 0.0366 Inf   0.423  0.9987
# pronoun - gap island      KE      -0.000382 0.0361 Inf  -0.011  1.0000
# pronoun - gap short       ME       0.025186 0.0319 Inf   0.789  0.9644
# pronoun - gap long        ME       0.028438 0.0324 Inf   0.877  0.9414
# pronoun - gap island      ME      -0.038343 0.0322 Inf  -1.191  0.7930

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_logrt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_logrt_proficiency_posthoc.csv') 
dt <- read_rds('tables/orc_ajt_logrt_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

tic() 
model2 <- lmer(log(rt) ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency + environment | participant) +
                 (1 + pairing | item),
               data = filter(md, study == 'SRC'), 
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 8 sec
# singular fit

check_convergence(model2) # okay
check_singularity(model2) # warning, but random-effects table looks okay

# set as final model ...

model <- model2

# save ...

model %>% write_rds('models/src_ajt_logrt_proficiency_model.rds')
model <- read_rds('models/src_ajt_logrt_proficiency_model.rds')

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_logrt_proficiency_model.rds') 
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_logrt_proficiency_model.csv') 
dt <- read_rds('tables/src_ajt_logrt_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                             |Est.  |$SE$ |df      |$t$    |$p$    |
# |:-----------------------------------------------------|:-----|:----|:-------|:------|:------|
# |(Intercept)                                           |8.82  |0.03 |133.44  |311.49 |< .001 | ***
# |Dependency[RP]                                        |0.05  |0.01 |140.28  |3.76   |< .001 | ***
# |Environment[Long]                                     |0.07  |0.02 |237.77  |4.01   |< .001 | ***
# |Environment[Island]                                   |0.13  |0.02 |133.58  |6.20   |< .001 | ***
# |Pairing[ME]                                           |0.18  |0.06 |134.52  |3.22   |0.002  | ***
# |Proficiency                                           |0.03  |0.03 |131.32  |1.12   |0.263  |
# |Dependency[RP]:Environment[Long]                      |−0.15 |0.03 |3245.02 |−4.84  |< .001 | ***
# |Dependency[RP]:Environment[Island]                    |−0.21 |0.03 |3353.35 |−6.60  |< .001 | ***
# |Dependency[RP]:Pairing[ME]                            |0.02  |0.03 |140.28  |0.73   |0.470  |
# |Environment[Long]:Pairing[ME]                         |−0.06 |0.03 |237.74  |−1.66  |0.098  |
# |Environment[Island]:Pairing[ME]                       |0.00  |0.04 |133.57  |−0.02  |0.984  |
# |Dependency[RP]:Proficiency                            |0.03  |0.01 |138.04  |1.89   |0.061  |
# |Environment[Long]:Proficiency                         |−0.01 |0.02 |230.25  |−0.35  |0.729  |
# |Environment[Island]:Proficiency                       |0.02  |0.02 |133.73  |1.03   |0.306  |
# |Pairing[ME]:Proficiency                               |−0.01 |0.06 |131.32  |−0.17  |0.865  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]               |0.01  |0.06 |3246.60 |0.19   |0.849  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME]             |0.04  |0.06 |3353.49 |0.59   |0.553  |
# |Dependency[RP]:Envir.[Long]:Proficiency               |0.05  |0.03 |3328.51 |1.55   |0.122  |
# |Dependency[RP]:Envir.[Island]:Proficiency             |0.02  |0.03 |3349.58 |0.47   |0.639  |
# |Dependency[RP]:Pairing[ME]:Proficiency                |0.01  |0.03 |138.03  |0.18   |0.854  |
# |Environment[Long]:Pairing[ME]:Proficiency             |−0.02 |0.03 |230.25  |−0.46  |0.645  |
# |Environment[Island]:Pairing[ME]:Proficiency           |0.00  |0.04 |133.73  |−0.07  |0.947  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]:Proficiency   |−0.05 |0.06 |3328.84 |−0.79  |0.430  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME]:Proficiency |−0.04 |0.06 |3349.51 |−0.61  |0.541  |

# joint tests ...

joint_tests(model)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 Inf  14.240  0.0002 ***
# environment                                  2 Inf  19.562  <.0001 ***
# pairing                                      1 Inf  10.359  0.0013 **
# proficiency                                  1 Inf   1.262  0.2612
# dependency:environment                       2 Inf  23.369  <.0001 ***
# dependency:pairing                           1 Inf   0.528  0.4676
# dependency:proficiency                       1 Inf   3.556  0.0593 .
# environment:pairing                          2 Inf   1.940  0.1436
# environment:proficiency                      2 Inf   1.109  0.3299
# pairing:proficiency                          1 Inf   0.029  0.8647
# dependency:environment:pairing               2 Inf   0.181  0.8342
# dependency:environment:proficiency           2 Inf   1.252  0.2860
# dependency:pairing:proficiency               1 Inf   0.034  0.8537
# environment:pairing:proficiency              2 Inf   0.130  0.8779
# dependency:environment:pairing:proficiency   2 Inf   0.344  0.7092

# post-hoc tests ...

pairwise <- model %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate     SE  df z.ratio p.value
# pronoun - gap short       KE       -0.0118 0.0341 Inf  -0.346  0.9996
# pronoun - gap long        KE        0.0625 0.0342 Inf   1.829  0.3411
# pronoun - gap island      KE        0.0229 0.0347 Inf   0.660  0.9858
# pronoun - gap short       ME        0.0233 0.0318 Inf   0.735  0.9757
# pronoun - gap long        ME        0.0474 0.0318 Inf   1.490  0.5839
# pronoun - gap island      ME        0.0188 0.0323 Inf   0.583  0.9927

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_logrt_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_logrt_proficiency_posthoc.csv') 
dt <- read_rds('tables/src_ajt_logrt_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: acceptance variance ####
# ---------------------------------------------------------------------------- #

# filter ...

crit <- ajt_clean %>% 
  filter(!condition %in% c('grammatical', 'ungrammatical'))

# calculate by-participant variances ...

ajt_var <- crit %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(', pairing, '; *n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, participant, task, dependency, environment, proficiency) %>%
  summarise(var = var(acceptance)) %>%
  ungroup()

# summarize ...

plot <- ajt_var %>%
  group_by(study, pairing, pairing_label, task, dependency, environment) %>%
  summarise(mean = mean(var, na.rm = T),
            lcl = my.lcl(var),
            ucl = my.ucl(var)) %>%
  ungroup() %>%
  filter(dependency != 'NA') %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = case_when(environment == 'Island' & pairing %in% c('KK', 'MM') ~ "‘Island’",
                                 TRUE ~ environment)) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island", "‘Island’")) %>%
  arrange(study, pairing, dependency, environment)

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, 
      shape=dependency),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                width = 0.2, 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name='Environment'),
  scale_y_continuous(name='Variance', breaks = c(0, .1, .2)),
  coord_cartesian(ylim = c(0, 0.2)),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
)

p1 + s 
ggsave("plots/orc_ajt_var.png", width=6.5, height=2.25, dpi=600)

p2 + s 
ggsave("plots/src_ajt_var.png", width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot %>%
  select(study, pairing, dependency, environment, mean, lcl, ucl) %>%
  arrange(study, pairing, dependency, environment) %>%
  my_table() 

dt %>% write_rds('tables/ajt_var.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/ajt_var.csv')
dt <- read_rds('tables/ajt_var.rds')

# ---------------------------------------------------------------------------- #
# model: acceptance variance ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- ajt_var %>%
  select(study, pairing, participant, dependency, environment, var, proficiency) %>%
  mutate_at(c('study', 'pairing', 'participant', 'dependency', 'environment'), factor) %>%
  mutate(var = as.numeric(var)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$pairing) <- contr.treatment(5) - matrix(rep(1/5, 20), ncol = 4)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model2 <- lmer(var ~ dependency * environment * pairing + 
                 (1 + dependency * environment | participant), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 12 sec
# singular fit

check_convergence(model2) # okay
check_singularity(model2) # warning, but random-effects table looks okay

anova(model1, model2, refit = FALSE) # model2 is not better (p = .944)

# set as final model

model <- model2

# save ...

model %>% write_rds('models/orc_ajt_var_model.rds')
model <- read_rds('models/orc_ajt_var_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_ajt_var_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_var_model.csv')
dt <- read_rds('tables/orc_ajt_var_model.rds')
csv <- read_csv('tables/orc_ajt_var_model.csv')

dt %>% rename("Est." = 2) %>% kbl(format = 'markdown')

# |Parameter                                 |Est.  |$SE$ |$t$   |$p$    |
# |:-----------------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                               |0.10  |0.00 |27.82 |< .001 | ***
# |Dependency[RP]                            |0.01  |0.01 |1.09  |0.278  |
# |Environment[Long]                         |0.02  |0.01 |3.10  |0.002  | **
# |Environment[Island]                       |0.03  |0.01 |5.73  |< .001 | ***
# |Pairing[KE]                               |0.03  |0.01 |2.98  |0.003  | **
# |Pairing[ME]                               |0.05  |0.01 |4.81  |< .001 | ***
# |Pairing[KK]                               |0.03  |0.01 |2.94  |0.004  | **
# |Pairing[MM]                               |0.06  |0.01 |5.95  |< .001 | ***
# |Dependency[RP]:Environment[Long]          |−0.02 |0.01 |−2.06 |0.040  | *
# |Dependency[RP]:Environment[Island]        |−0.02 |0.01 |−1.91 |0.056  | .
# |Dependency[RP]:Pairing[KE]                |−0.02 |0.02 |−1.00 |0.320  |
# |Dependency[RP]:Pairing[ME]                |−0.01 |0.02 |−0.50 |0.616  |
# |Dependency[RP]:Pairing[KK]                |0.01  |0.02 |0.37  |0.710  |
# |Dependency[RP]:Pairing[MM]                |−0.02 |0.02 |−1.29 |0.198  |
# |Environment[Long]:Pairing[KE]             |−0.01 |0.02 |−0.42 |0.673  |
# |Environment[Island]:Pairing[KE]           |−0.01 |0.02 |−0.67 |0.500  |
# |Environment[Long]:Pairing[ME]             |−0.01 |0.02 |−0.47 |0.635  |
# |Environment[Island]:Pairing[ME]           |−0.01 |0.02 |−0.53 |0.599  |
# |Environment[Long]:Pairing[KK]             |0.03  |0.02 |1.44  |0.151  |
# |Environment[Island]:Pairing[KK]           |0.04  |0.02 |2.45  |0.015  | *
# |Environment[Long]:Pairing[MM]             |−0.03 |0.02 |−1.68 |0.092  | .
# |Environment[Island]:Pairing[MM]           |−0.05 |0.02 |−2.62 |0.009  | **
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |0.02  |0.04 |0.61  |0.544  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |0.07  |0.04 |1.86  |0.063  | .
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |0.03  |0.04 |0.94  |0.346  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |0.04  |0.04 |1.10  |0.272  |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |0.00  |0.04 |−0.03 |0.979  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−0.02 |0.04 |−0.45 |0.652  |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |0.08  |0.04 |2.26  |0.024  | *
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |0.10  |0.04 |2.85  |0.004  | **

# joint tests ...

joint_tests(model)

# model term                     df1    df2 F.ratio p.value
# dependency                       1 195.99   1.184  0.2779
# environment                      2 195.70  16.429  <.0001 ***
# pairing                          4 357.26   9.073  <.0001 ***
# dependency:environment           2 195.69   2.654  0.0729 .
# dependency:pairing               4 482.77   0.934  0.4440
# environment:pairing              8 570.29   3.021  0.0025 **
# dependency:environment:pairing   8 566.07   1.921  0.0546 .

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing  estimate     SE  df t.ratio p.value
# pronoun - gap short       EE       0.050562 0.0174 529   2.899  0.0567 .
# pronoun - gap long        EE       0.000000 0.0176 519   0.000  1.0000
# pronoun - gap island      EE      -0.010112 0.0170 561  -0.595  1.0000
# pronoun - gap short       KE       0.004774 0.0200 531   0.239  1.0000
# pronoun - gap long        KE      -0.023880 0.0201 521  -1.186  0.9814
# pronoun - gap island      KE       0.011802 0.0194 561   0.607  1.0000
# pronoun - gap short       ME       0.018860 0.0189 529   0.999  0.9966
# pronoun - gap long        ME       0.001316 0.0190 519   0.069  1.0000
# pronoun - gap island      ME      -0.003070 0.0184 561  -0.167  1.0000
# pronoun - gap short       KK       0.062301 0.0200 531   3.122  0.0279 *
# pronoun - gap long        KK       0.010789 0.0201 521   0.536  1.0000
# pronoun - gap island      KK      -0.014788 0.0194 561  -0.761  0.9998
# pronoun - gap short       MM      -0.029458 0.0190 530  -1.550  0.8532
# pronoun - gap long        MM      -0.000799 0.0192 520  -0.042  1.0000
# pronoun - gap island      MM       0.010630 0.0185 561   0.574  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_var_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_var_posthoc.csv')
dt <- read_rds('tables/orc_ajt_var_posthoc.rds')
csv <- read_csv('tables/orc_ajt_var_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model2 <- lmer(var ~ dependency * environment * pairing + 
                 (1 + dependency * environment | participant), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 14 sec
# singular fit

check_convergence(model2) # okay
check_singularity(model2) # warning, but random-effects table looks okay

anova(model1, model2, refit = FALSE) # model2 is not better

# set as final model

model <- model2

# save ...

model %>% write_rds('models/src_ajt_var_model.rds')
model <- read_rds('models/src_ajt_var_model.rds')

# summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_ajt_var_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_var_model.csv')
dt <- read_rds('tables/src_ajt_var_model.rds')
csv <- read_csv('tables/src_ajt_var_model.csv')

dt %>% rename("Est." = 2) %>% kbl(format = 'markdown')

# |Parameter                                 |Est.  |$SE$ |$t$   |$p$    |
# |:-----------------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                               |0.09  |0.00 |24.07 |< .001 | ***
# |Dependency[RP]                            |0.03  |0.00 |5.28  |< .001 | ***
# |Environment[Long]                         |0.04  |0.01 |6.74  |< .001 | ***
# |Environment[Island]                       |0.05  |0.01 |7.35  |< .001 | ***
# |Pairing[KE]                               |0.04  |0.01 |4.14  |< .001 | ***
# |Pairing[ME]                               |0.07  |0.01 |6.33  |< .001 | ***
# |Pairing[KK]                               |0.04  |0.01 |3.48  |< .001 | ***
# |Pairing[MM]                               |0.07  |0.01 |7.11  |< .001 | ***
# |Dependency[RP]:Environment[Long]          |−0.05 |0.01 |−4.16 |< .001 | ***
# |Dependency[RP]:Environment[Island]        |−0.04 |0.01 |−3.17 |0.002  | **
# |Dependency[RP]:Pairing[KE]                |0.01  |0.02 |0.44  |0.657  |
# |Dependency[RP]:Pairing[ME]                |−0.01 |0.02 |−0.59 |0.558  |
# |Dependency[RP]:Pairing[KK]                |0.01  |0.02 |0.70  |0.485  |
# |Dependency[RP]:Pairing[MM]                |0.02  |0.02 |1.07  |0.283  |
# |Environment[Long]:Pairing[KE]             |0.02  |0.02 |0.79  |0.428  |
# |Environment[Island]:Pairing[KE]           |0.02  |0.02 |1.10  |0.270  |
# |Environment[Long]:Pairing[ME]             |−0.02 |0.02 |−1.22 |0.221  |
# |Environment[Island]:Pairing[ME]           |−0.03 |0.02 |−1.43 |0.153  |
# |Environment[Long]:Pairing[KK]             |0.02  |0.02 |0.87  |0.385  |
# |Environment[Island]:Pairing[KK]           |0.03  |0.02 |1.67  |0.096  | .
# |Environment[Long]:Pairing[MM]             |−0.01 |0.02 |−0.36 |0.718  |
# |Environment[Island]:Pairing[MM]           |0.01  |0.02 |0.58  |0.559  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |0.00  |0.04 |−0.12 |0.903  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |−0.02 |0.04 |−0.56 |0.573  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |−0.07 |0.04 |−1.81 |0.071  | .
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |−0.07 |0.04 |−1.80 |0.073  | .
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |0.00  |0.04 |0.11  |0.910  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−0.09 |0.04 |−2.23 |0.026  | *
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |−0.09 |0.04 |−2.27 |0.023  | *
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |−0.13 |0.04 |−3.34 |< .001 | ***

# joint tests ...

joint_tests(model)

# model term                     df1    df2 F.ratio p.value
# dependency                       1 189.42  27.913  <.0001 ***
# environment                      2 187.72  32.774  <.0001 ***
# pairing                          4 303.24  13.375  <.0001 ***
# dependency:environment           2 188.98   9.282  0.0001 ***
# dependency:pairing               4 482.62   0.964  0.4268
# environment:pairing              8 445.33   1.711  0.0935 .
# dependency:environment:pairing   8 477.28   2.783  0.0051 **

# post-hoc tests ...

pairwise <- model %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate     SE  df t.ratio p.value
# pronoun - gap short       EE       0.02038 0.0197 500   1.036  0.9950
# pronoun - gap long        EE       0.00087 0.0198 490   0.044  1.0000
# pronoun - gap island      EE       0.04242 0.0201 476   2.112  0.4132
# pronoun - gap short       KE       0.03636 0.0187 497   1.948  0.5472
# pronoun - gap long        KE       0.01212 0.0188 487   0.644  1.0000
# pronoun - gap island      KE       0.03636 0.0190 476   1.915  0.5757
# pronoun - gap short       ME       0.05699 0.0179 498   3.189  0.0225 *
# pronoun - gap long        ME      -0.03098 0.0180 488  -1.719  0.7379
# pronoun - gap island      ME       0.01028 0.0182 478   0.565  1.0000
# pronoun - gap short       KK       0.05909 0.0187 497   3.166  0.0243 *
# pronoun - gap long        KK       0.04394 0.0188 487   2.334  0.2600
# pronoun - gap island      KK      -0.00606 0.0190 476  -0.319  1.0000
# pronoun - gap short       MM       0.10822 0.0177 497   6.098  <.0001 ***
# pronoun - gap long        MM       0.00274 0.0179 487   0.153  1.0000
# pronoun - gap island      MM       0.00274 0.0181 476   0.152  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_var_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_var_posthoc.csv')
dt <- read_rds('tables/src_ajt_var_posthoc.rds')
csv <- read_csv('tables/src_ajt_var_posthoc.csv')

# plot ...

plot(emmeans(model, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: acceptance variance proficiency ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- ajt_var %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(is.na(proficiency) == FALSE) %>%
  select(study, pairing, participant, dependency, environment, var, proficiency) %>%
  mutate_at(c('study', 'pairing', 'participant', 'dependency', 'environment'), factor) %>%
  mutate(var = as.numeric(var)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$pairing) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(var ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 2 sec
# no warnings

check_convergence(model1) # okay
check_singularity(model1) # okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/orc_ajt_var_proficiency_model.rds')
model <- read_rds('models/orc_ajt_var_proficiency_model.rds')

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_var_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_var_proficiency_model.csv')
dt <- read_rds('tables/orc_ajt_var_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                             |Est.  |$SE$ |df     |$t$   |$p$    |
# |:-----------------------------------------------------|:-----|:----|:------|:-----|:------|
# |(Intercept)                                           |0.11  |0.01 |140.00 |19.13 |< .001 | ***
# |Dependency[RP]                                        |0.00  |0.01 |140.00 |0.20  |0.838  |
# |Environment[Long]                                     |0.01  |0.01 |140.00 |1.57  |0.118  |
# |Environment[Island]                                   |0.03  |0.01 |140.00 |2.75  |0.007  | **
# |Pairing[ME]                                           |0.01  |0.01 |140.00 |1.29  |0.198  |
# |Proficiency                                           |−0.03 |0.01 |140.00 |−5.12 |< .001 | ***
# |Dependency[RP]:Environment[Long]                      |−0.02 |0.02 |280.00 |−1.37 |0.173  |
# |Dependency[RP]:Environment[Island]                    |−0.01 |0.02 |280.00 |−0.45 |0.650  |
# |Dependency[RP]:Pairing[ME]                            |0.01  |0.02 |140.00 |0.47  |0.642  |
# |Environment[Long]:Pairing[ME]                         |0.00  |0.02 |140.00 |−0.05 |0.963  |
# |Environment[Island]:Pairing[ME]                       |0.00  |0.02 |140.00 |0.18  |0.860  |
# |Dependency[RP]:Proficiency                            |0.01  |0.01 |140.00 |1.20  |0.232  |
# |Environment[Long]:Proficiency                         |0.00  |0.01 |140.00 |−0.23 |0.816  |
# |Environment[Island]:Proficiency                       |0.02  |0.01 |140.00 |1.50  |0.136  |
# |Pairing[ME]:Proficiency                               |0.02  |0.01 |140.00 |1.42  |0.156  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]               |0.01  |0.03 |280.00 |0.34  |0.732  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME]             |−0.03 |0.03 |280.00 |−0.85 |0.396  |
# |Dependency[RP]:Envir.[Long]:Proficiency               |0.01  |0.02 |280.00 |0.56  |0.578  |
# |Dependency[RP]:Envir.[Island]:Proficiency             |0.00  |0.02 |280.00 |0.15  |0.882  |
# |Dependency[RP]:Pairing[ME]:Proficiency                |0.01  |0.02 |140.00 |0.52  |0.603  |
# |Environment[Long]:Pairing[ME]:Proficiency             |0.01  |0.02 |140.00 |0.31  |0.755  |
# |Environment[Island]:Pairing[ME]:Proficiency           |−0.02 |0.02 |140.00 |−1.04 |0.301  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]:Proficiency   |0.02  |0.03 |280.00 |0.70  |0.486  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME]:Proficiency |0.00  |0.03 |280.00 |0.14  |0.888  |

# joint tests ...

joint_tests(model)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 140   0.042  0.8382
# environment                                  2 140   3.779  0.0252 *
# pairing                                      1 140   1.671  0.1983
# proficiency                                  1 140  26.189  <.0001 ***
# dependency:environment                       2 280   0.968  0.3811
# dependency:pairing                           1 140   0.217  0.6422
# dependency:proficiency                       1 140   1.444  0.2315
# environment:pairing                          2 140   0.030  0.9700
# environment:proficiency                      2 140   1.939  0.1477
# pairing:proficiency                          1 140   2.030  0.1564
# dependency:environment:pairing               2 280   0.754  0.4713
# dependency:environment:proficiency           2 280   0.167  0.8467
# dependency:pairing:proficiency               1 140   0.272  0.6030
# environment:pairing:proficiency              2 140   1.105  0.3339
# dependency:environment:pairing:proficiency   2 280   0.272  0.7619

# post-hoc tests ...

pairwise <- model %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate     SE  df t.ratio p.value
# pronoun - gap short       KE       0.00692 0.0200 391   0.347  0.9996
# pronoun - gap long        KE       0.00454 0.0200 391   0.227  1.0000
# pronoun - gap island      KE       0.00703 0.0200 391   0.352  0.9995
# pronoun - gap short       ME       0.00690 0.0176 391   0.391  0.9991
# pronoun - gap long        ME       0.02811 0.0176 391   1.595  0.4979
# pronoun - gap island      ME       0.01179 0.0176 391   0.669  0.9837

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_var_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_var_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_ajt_var_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model1 <- lmer(var ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 2 sec
# singular fit

check_convergence(model1) # okay
check_singularity(model1) # warning, but random-effects table looks okay

# set as final model ...

model <- model1

# save ...

model %>% write_rds('models/src_ajt_var_proficiency_model.rds')
model <- read_rds('models/src_ajt_var_proficiency_model.rds')

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_var_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_var_proficiency_model.csv')
dt <- read_rds('tables/src_ajt_var_proficiency_model.rds')

dt %>% rename("Est." = 2) %>% kable(format = "markdown")

# |Parameter                                             |Est.  |$SE$ |df     |$t$   |$p$    |
# |:-----------------------------------------------------|:-----|:----|:------|:-----|:------|
# |(Intercept)                                           |0.10  |0.01 |134.10 |19.08 |< .001 | ***
# |Dependency[RP]                                        |0.02  |0.01 |157.47 |2.43  |0.016  | *
# |Environment[Long]                                     |0.04  |0.01 |259.70 |3.82  |< .001 | ***
# |Environment[Island]                                   |0.03  |0.01 |147.02 |3.50  |< .001 | ***
# |Pairing[ME]                                           |0.02  |0.01 |134.10 |1.41  |0.162  |
# |Proficiency                                           |−0.01 |0.01 |134.10 |−2.47 |0.015  | *
# |Dependency[RP]:Environment[Long]                      |−0.07 |0.02 |399.00 |−3.68 |< .001 | ***
# |Dependency[RP]:Environment[Island]                    |−0.03 |0.02 |399.00 |−1.60 |0.110  |
# |Dependency[RP]:Pairing[ME]                            |−0.01 |0.02 |157.47 |−0.86 |0.391  |
# |Environment[Long]:Pairing[ME]                         |−0.04 |0.02 |259.70 |−1.92 |0.056  | .
# |Environment[Island]:Pairing[ME]                       |−0.05 |0.02 |147.02 |−2.55 |0.012  | *
# |Dependency[RP]:Proficiency                            |0.00  |0.01 |157.47 |0.13  |0.898  |
# |Environment[Long]:Proficiency                         |0.00  |0.01 |259.70 |0.02  |0.987  |
# |Environment[Island]:Proficiency                       |0.00  |0.01 |147.02 |−0.28 |0.779  |
# |Pairing[ME]:Proficiency                               |0.03  |0.01 |134.10 |2.43  |0.017  | *
# |Dependency[RP]:Envir.[Long]:Pairing[ME]               |−0.07 |0.04 |399.00 |−1.82 |0.069  | .
# |Dependency[RP]:Envir.[Island]:Pairing[ME]             |−0.04 |0.04 |399.00 |−1.17 |0.243  |
# |Dependency[RP]:Envir.[Long]:Proficiency               |−0.01 |0.02 |399.00 |−0.49 |0.624  |
# |Dependency[RP]:Envir.[Island]:Proficiency             |0.01  |0.02 |399.00 |0.39  |0.697  |
# |Dependency[RP]:Pairing[ME]:Proficiency                |0.00  |0.02 |157.47 |−0.27 |0.790  |
# |Environment[Long]:Pairing[ME]:Proficiency             |−0.01 |0.02 |259.70 |−0.61 |0.544  |
# |Environment[Island]:Pairing[ME]:Proficiency           |−0.03 |0.02 |147.02 |−1.31 |0.193  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]:Proficiency   |−0.13 |0.04 |399.00 |−3.39 |< .001 | ***
# |Dependency[RP]:Envir.[Island]:Pairing[ME]:Proficiency |−0.07 |0.04 |399.00 |−1.83 |0.069  | .

# joint tests ...

joint_tests(model)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 133   5.885  0.0166 *
# environment                                  2 133   8.942  0.0002 ***
# pairing                                      1 133   1.979  0.1618
# proficiency                                  1 133   6.102  0.0148 *
# dependency:environment                       2 266   6.813  0.0013 **
# dependency:pairing                           1 133   0.739  0.3914
# dependency:proficiency                       1 133   0.017  0.8979
# environment:pairing                          2 133   3.519  0.0324 *
# environment:proficiency                      2 133   0.056  0.9453
# pairing:proficiency                          1 133   5.887  0.0166 *
# dependency:environment:pairing               2 266   1.706  0.1835
# dependency:environment:proficiency           2 266   0.389  0.6779
# dependency:pairing:proficiency               1 133   0.071  0.7902
# environment:pairing:proficiency              2 133   0.857  0.4267
# dependency:environment:pairing:proficiency   2 266   5.764  0.0035 **

# post-hoc tests ...

pairwise <- model %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate     SE  df t.ratio p.value
# pronoun - gap short       KE      -0.02888 0.0203 396  -1.426  0.6325
# pronoun - gap long        KE       0.02593 0.0203 396   1.280  0.7378
# pronoun - gap island      KE       0.01298 0.0203 396   0.641  0.9878
# pronoun - gap short       ME       0.03233 0.0188 396   1.724  0.4131
# pronoun - gap long        ME      -0.04103 0.0188 396  -2.187  0.1627
# pronoun - gap island      ME       0.00519 0.0188 396   0.277  0.9999

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_var_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_var_proficiency_posthoc.csv')
dt <- read_rds('tables/src_ajt_var_proficiency_posthoc.rds')

# plot ... 

emmip(model, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: acceptance variance proficiency ####
# ---------------------------------------------------------------------------- #

plot <- md %>%
  group_by(study, pairing) %>%
  mutate(group = case_when(pairing == 'KE' ~ 'KLE',
                           pairing == 'ME' ~ 'MLE')) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, '<br>(*n* = ', pairing_size, ')')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = var, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='Variance'),
  coord_cartesian(ylim = c(0, .3)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ajt_var_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s
ggsave('plots/src_ajt_var_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: iqr ####
# ---------------------------------------------------------------------------- #

# filter ...

crit <- ajt_clean %>% 
  filter(!condition %in% c('grammatical', 'ungrammatical'))

# calculate by-participant iqr values ...

ajt_iqr <- crit %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, participant, task, dependency, environment, proficiency) %>%
  summarise(iqr = IQR(response)) %>%
  ungroup()

# summarize ...
  
plot <- ajt_iqr %>%
  group_by(study, pairing, pairing_label, task, dependency, environment) %>%
  summarise(mean = mean(iqr, na.rm = T),
            lcl = my.lcl(iqr),
            ucl = my.ucl(iqr)) %>%
  ungroup() %>%
  filter(dependency != 'NA') %>%
  mutate(environment = case_when(environment == 'short' ~ 'Short',
                                 environment == 'long' ~ 'Long',
                                 environment == 'island' ~ 'Island')) %>%
  mutate(environment = case_when(environment == 'Island' & pairing %in% c('KK', 'MM') ~ "‘Island’",
                                 TRUE ~ environment)) %>%
  mutate(environment = fct_relevel(environment, 
                                   "Short", "Long", "Island", "‘Island’")) %>%
  arrange(study, pairing, dependency, environment)

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))

p2 <- ggplot(data=filter(plot, study == 'SRC'))

design <- '
 ABD
 #CE
'

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, 
      shape=dependency),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                width = 0.2, 
                linetype = 1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_x_discrete(name='Environment'),
  scale_y_continuous(name='IQR'),
  coord_cartesian(ylim = c(0, 1.7)),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency', 
                        values=c('41', 'solid'), 
                        labels=c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'bottom',
        legend.margin = margin(t = -.4, unit = 'cm'),
        panel.spacing = unit(.1, "lines"),
        plot.margin = margin(0, 0, 0, 0)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x')
  # theme(
  #   text = element_text(size = 12),
  #   plot.margin = margin(2, 2, 2, 2),
  #   legend.position = c(.15, .22),
  #   legend.margin = margin(-10, 0, 0, 0),
  #   axis.title.y = element_text(hjust = .88),
  #   strip.text = element_markdown()),
  # facet_manual(~pairing_label, design = design, axes = 'all', 
  #              remove_labels = FALSE, scales = 'free')
)

p1 + s 
ggsave("plots/orc_ajt_iqr.png", width=6.5, height=2.25, dpi=600)

p2 + s
ggsave("plots/src_ajt_iqr.png", width=6.5, height=2.25, dpi=600)

# table ...

dt <- plot %>%
  select(study, pairing, dependency, environment, mean, lcl, ucl) %>%
  arrange(study, pairing, dependency, environment) %>%
  my_table() 

dt %>% write_rds('tables/ajt_iqr.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/ajt_iqr.csv')
dt <- read_rds('tables/ajt_iqr.rds')

# ---------------------------------------------------------------------------- #
# model: iqr ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- ajt_iqr %>%
  select(study, pairing, participant, dependency, environment, iqr) %>%
  mutate_at(c('study', 'pairing', 'participant', 'dependency', 'environment'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$pairing) <- contr.treatment(5) - matrix(rep(1/5, 20), ncol = 4) 

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit ...

tic()
model1 <- lmer(iqr ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# fast

check_convergence(model1) # okay
check_singularity(model1) # warning (but random effects table looks okay)

# save ...

model1 %>% write_rds('models/orc_ajt_iqr_model.rds')
model1 <- read_rds('models/orc_ajt_iqr_model.rds')

# summary ...

summary(model1)

# Parameter                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                          0.81606    0.03542  221.29621  23.040  < 2e-16 ***
# dependency2                          0.05376    0.05586  214.99096   0.962  0.33695    
# environment2                         0.11480    0.06211  649.11407   1.848  0.06502 .  
# environment3                         0.28898    0.06339  498.76292   4.559 6.49e-06 ***
# pairing2                             0.24816    0.09657  425.95565   2.570  0.01051 *  
# pairing3                             0.40702    0.09370  422.99561   4.344 1.75e-05 ***
# pairing4                             0.16374    0.09657  425.95565   1.696  0.09069 .  
# pairing5                             0.46642    0.09399  425.56675   4.963 1.01e-06 ***
# dependency2:environment2            -0.01825    0.12152 1630.67544  -0.150  0.88066    
# dependency2:environment3            -0.03005    0.12152 1630.67544  -0.247  0.80470    
# dependency2:pairing2                -0.07933    0.16512  551.30869  -0.480  0.63109    
# dependency2:pairing3                -0.08587    0.16013  549.55591  -0.536  0.59200    
# dependency2:pairing4                 0.18959    0.16512  551.30869   1.148  0.25140    
# dependency2:pairing5                -0.13012    0.16069  551.07923  -0.810  0.41843    
# environment2:pairing2               -0.05149    0.19116 1206.13867  -0.269  0.78769    
# environment3:pairing2               -0.07734    0.19337 1035.58912  -0.400  0.68927    
# environment2:pairing3               -0.08560    0.18537 1205.31084  -0.462  0.64432    
# environment3:pairing3               -0.10441    0.18752 1033.97889  -0.557  0.57778    
# environment2:pairing4               -0.03238    0.19116 1206.13867  -0.169  0.86554    
# environment3:pairing4                0.24583    0.19337 1035.58912   1.271  0.20392    
# environment2:pairing5               -0.40823    0.18603 1206.03052  -2.194  0.02840 *  
# environment3:pairing5               -0.37479    0.18819 1035.37869  -1.992  0.04668 *  
# dependency2:environment2:pairing2    0.57506    0.37769 1630.67544   1.523  0.12806    
# dependency2:environment3:pairing2    0.67044    0.37769 1630.67544   1.775  0.07607 .  
# dependency2:environment2:pairing3    0.98315    0.36625 1630.67544   2.684  0.00734 ** 
# dependency2:environment3:pairing3    0.76796    0.36625 1630.67544   2.097  0.03616 *  
# dependency2:environment2:pairing4    0.23315    0.37769 1630.67544   0.617  0.53713    
# dependency2:environment3:pairing4    0.36897    0.37769 1630.67544   0.977  0.32876    
# dependency2:environment2:pairing5    0.53315    0.36756 1630.67544   1.450  0.14711    
# dependency2:environment3:pairing5    0.40191    0.36756 1630.67544   1.093  0.27436   
  
# table ...

dt <- summary(model1)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/orc_ajt_iqr_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_iqr_model.csv')
dt <- read_rds('tables/orc_ajt_iqr_model.rds')
csv <- read_csv('tables/orc_ajt_iqr_model.csv')

dt %>% kbl(format = 'markdown')

# |Parameter                                 |Estimate~($\hat{\beta}$) |$SE$ |$t$   |$p$    |
# |:-----------------------------------------|:------------------------|:----|:-----|:------|
# |(Intercept)                               |0.82                     |0.04 |23.04 |< .001 |
# |Dependency[RP]                            |0.05                     |0.06 |0.96  |0.337  |
# |Environment[Long]                         |0.11                     |0.06 |1.85  |0.065  |
# |Environment[Island]                       |0.29                     |0.06 |4.56  |< .001 |
# |Pairing[KE]                               |0.25                     |0.10 |2.57  |0.011  |
# |Pairing[ME]                               |0.41                     |0.09 |4.34  |< .001 |
# |Pairing[KK]                               |0.16                     |0.10 |1.70  |0.091  |
# |Pairing[MM]                               |0.47                     |0.09 |4.96  |< .001 |
# |Dependency[RP]:Environment[Long]          |−0.02                    |0.12 |−0.15 |0.881  |
# |Dependency[RP]:Environment[Island]        |−0.03                    |0.12 |−0.25 |0.805  |
# |Dependency[RP]:Pairing[KE]                |−0.08                    |0.17 |−0.48 |0.631  |
# |Dependency[RP]:Pairing[ME]                |−0.09                    |0.16 |−0.54 |0.592  |
# |Dependency[RP]:Pairing[KK]                |0.19                     |0.17 |1.15  |0.251  |
# |Dependency[RP]:Pairing[MM]                |−0.13                    |0.16 |−0.81 |0.418  |
# |Environment[Long]:Pairing[KE]             |−0.05                    |0.19 |−0.27 |0.788  |
# |Environment[Island]:Pairing[KE]           |−0.08                    |0.19 |−0.40 |0.689  |
# |Environment[Long]:Pairing[ME]             |−0.09                    |0.19 |−0.46 |0.644  |
# |Environment[Island]:Pairing[ME]           |−0.10                    |0.19 |−0.56 |0.578  |
# |Environment[Long]:Pairing[KK]             |−0.03                    |0.19 |−0.17 |0.866  |
# |Environment[Island]:Pairing[KK]           |0.25                     |0.19 |1.27  |0.204  |
# |Environment[Long]:Pairing[MM]             |−0.41                    |0.19 |−2.19 |0.028  |
# |Environment[Island]:Pairing[MM]           |−0.37                    |0.19 |−1.99 |0.047  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |0.58                     |0.38 |1.52  |0.128  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |0.67                     |0.38 |1.78  |0.076  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |0.98                     |0.37 |2.68  |0.007  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |0.77                     |0.37 |2.10  |0.036  |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |0.23                     |0.38 |0.62  |0.537  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |0.37                     |0.38 |0.98  |0.329  |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |0.53                     |0.37 |1.45  |0.147  |
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |0.40                     |0.37 |1.09  |0.274  |

# joint tests ...

joint_tests(model1)

# model term                     df1     df2 F.ratio p.value
# dependency                       1  195.81   0.926  0.3371
# environment                      2  195.61  10.390  0.0001 ***
# pairing                          4  386.18   7.006  <.0001 ***
# dependency:environment           2 1359.13   0.031  0.9694
# dependency:pairing               4  505.91   1.105  0.3533
# environment:pairing              8  562.53   1.614  0.1177
# dependency:environment:pairing   8 1359.13   1.215  0.2863

# post-hoc tests ...

pairwise <- model1 %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate    SE   df t.ratio p.value
# pronoun - gap short       EE       0.39326 0.180 1464   2.185  0.3564
# pronoun - gap long        EE      -0.08989 0.180 1464  -0.499  1.0000
# pronoun - gap island      EE      -0.07865 0.180 1464  -0.437  1.0000
# pronoun - gap short       KE      -0.10124 0.206 1466  -0.492  1.0000
# pronoun - gap long        KE      -0.00933 0.206 1466  -0.045  1.0000
# pronoun - gap island      KE       0.09729 0.206 1466   0.472  1.0000
# pronoun - gap short       ME      -0.27632 0.195 1464  -1.418  0.9205
# pronoun - gap long        ME       0.22368 0.195 1464   1.148  0.9866
# pronoun - gap island      ME       0.01974 0.195 1464   0.101  1.0000
# pronoun - gap short       KK       0.38214 0.206 1466   1.856  0.6254
# pronoun - gap long        KK       0.13214 0.206 1466   0.642  1.0000
# pronoun - gap island      KK       0.27920 0.206 1466   1.356  0.9435
# pronoun - gap short       MM      -0.04855 0.196 1465  -0.248  1.0000
# pronoun - gap long        MM       0.00145 0.196 1465   0.007  1.0000
# pronoun - gap island      MM      -0.11855 0.196 1465  -0.605  1.0000

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_iqr_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_iqr_posthoc.csv')
dt <- read_rds('tables/orc_ajt_iqr_posthoc.rds')
csv <- read_csv('tables/orc_ajt_iqr_posthoc.csv')

# plot ...

plot(emmeans(model1, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit ...

tic()
model1 <- lmer(iqr ~ dependency * environment * pairing + 
                 (1 + dependency + environment | participant), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model1)
toc()
beep()

# 2.34 sec

check_convergence(model1) # okay
check_singularity(model1) # warning (but random effects table looks okay)

# save ...

model1 %>% write_rds('models/src_ajt_iqr_model.rds')
model1 <- read_rds('models/src_ajt_iqr_model.rds')

# summary ...

summary(model1)

# (Intercept)                        7.730e-01  3.627e-02  2.073e+02  21.312  < 2e-16 ***
# dependency2                        1.862e-01  5.268e-02  2.796e+02   3.535 0.000477 ***
# environment2                       3.288e-01  6.373e-02  5.951e+02   5.159 3.38e-07 ***
# environment3                       4.180e-01  6.798e-02  3.701e+02   6.149 2.02e-09 ***
# pairing2                           3.965e-01  1.044e-01  3.803e+02   3.798 0.000169 ***
# pairing3                           4.904e-01  1.022e-01  3.824e+02   4.797 2.32e-06 ***
# pairing4                           2.626e-01  1.044e-01  3.803e+02   2.516 0.012276 *  
# pairing5                           6.196e-01  1.020e-01  3.804e+02   6.076 2.99e-09 ***
# dependency2:environment2          -5.991e-01  1.236e-01  1.651e+03  -4.848 1.37e-06 ***
# dependency2:environment3          -3.056e-01  1.236e-01  1.651e+03  -2.473 0.013499 *  
# dependency2:pairing2               1.817e-02  1.691e-01  7.115e+02   0.107 0.914474    
# dependency2:pairing3               8.751e-02  1.658e-01  7.124e+02   0.528 0.597682    
# dependency2:pairing4               3.465e-01  1.691e-01  7.115e+02   2.048 0.040885 *  
# dependency2:pairing5               3.426e-01  1.652e-01  7.117e+02   2.073 0.038511 *  
# environment2:pairing2              2.647e-01  2.056e-01  1.109e+03   1.288 0.198104    
# environment3:pairing2              2.701e-01  2.136e-01  7.844e+02   1.264 0.206482    
# environment2:pairing3             -7.052e-02  2.015e-01  1.110e+03  -0.350 0.726381    
# environment3:pairing3             -2.558e-01  2.093e-01  7.868e+02  -1.222 0.222112    
# environment2:pairing4              7.157e-03  2.056e-01  1.109e+03   0.035 0.972233    
# environment3:pairing4              5.798e-02  2.136e-01  7.844e+02   0.271 0.786158    
# environment2:pairing5             -1.124e-01  2.008e-01  1.109e+03  -0.560 0.575855    
# environment3:pairing5              3.494e-02  2.087e-01  7.845e+02   0.167 0.867092    
# dependency2:environment2:pairing2 -7.356e-02  4.040e-01  1.651e+03  -0.182 0.855521    
# dependency2:environment3:pairing2 -5.374e-02  4.040e-01  1.653e+03  -0.133 0.894194    
# dependency2:environment2:pairing3 -5.193e-01  3.959e-01  1.652e+03  -1.312 0.189814    
# dependency2:environment3:pairing3 -1.737e-01  3.959e-01  1.653e+03  -0.439 0.660948    
# dependency2:environment2:pairing4  7.795e-02  4.040e-01  1.651e+03   0.193 0.847009    
# dependency2:environment3:pairing4 -6.295e-01  4.040e-01  1.653e+03  -1.558 0.119396    
# dependency2:environment2:pairing5 -4.243e-01  3.947e-01  1.652e+03  -1.075 0.282443    
# dependency2:environment3:pairing5 -6.370e-01  3.947e-01  1.653e+03  -1.614 0.106764     

# table ...

dt <- summary(model1)$coefficients %>%
  as.data.frame() %>%
  select(-df) %>%
  rownames_to_column('Parameter') %>%
  my_table() 

dt %>% write_rds('tables/src_ajt_iqr_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_iqr_model.csv')
dt <- read_rds('tables/src_ajt_iqr_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                 |Estimate~($\hat{\beta}$) |$SE$ |$t$   |$p$    |
# |:-----------------------------------------|:------------------------|:----|:-----|:------|
# |(Intercept)                               |0.77                     |0.04 |21.31 |< .001 |
# |Dependency[RP]                            |0.19                     |0.05 |3.54  |< .001 |
# |Environment[Long]                         |0.33                     |0.06 |5.16  |< .001 |
# |Environment[Island]                       |0.42                     |0.07 |6.15  |< .001 |
# |Pairing[KE]                               |0.40                     |0.10 |3.80  |< .001 |
# |Pairing[ME]                               |0.49                     |0.10 |4.80  |< .001 |
# |Pairing[KK]                               |0.26                     |0.10 |2.52  |0.012  |
# |Pairing[MM]                               |0.62                     |0.10 |6.08  |< .001 |
# |Dependency[RP]:Environment[Long]          |−0.60                    |0.12 |−4.85 |< .001 |
# |Dependency[RP]:Environment[Island]        |−0.31                    |0.12 |−2.47 |0.013  |
# |Dependency[RP]:Pairing[KE]                |0.02                     |0.17 |0.11  |0.914  |
# |Dependency[RP]:Pairing[ME]                |0.09                     |0.17 |0.53  |0.598  |
# |Dependency[RP]:Pairing[KK]                |0.35                     |0.17 |2.05  |0.041  |
# |Dependency[RP]:Pairing[MM]                |0.34                     |0.17 |2.07  |0.039  |
# |Environment[Long]:Pairing[KE]             |0.26                     |0.21 |1.29  |0.198  |
# |Environment[Island]:Pairing[KE]           |0.27                     |0.21 |1.26  |0.206  |
# |Environment[Long]:Pairing[ME]             |−0.07                    |0.20 |−0.35 |0.726  |
# |Environment[Island]:Pairing[ME]           |−0.26                    |0.21 |−1.22 |0.222  |
# |Environment[Long]:Pairing[KK]             |0.01                     |0.21 |0.03  |0.972  |
# |Environment[Island]:Pairing[KK]           |0.06                     |0.21 |0.27  |0.786  |
# |Environment[Long]:Pairing[MM]             |−0.11                    |0.20 |−0.56 |0.576  |
# |Environment[Island]:Pairing[MM]           |0.03                     |0.21 |0.17  |0.867  |
# |Dependency[RP]:Envir.[Long]:Pairing[KE]   |−0.07                    |0.40 |−0.18 |0.856  |
# |Dependency[RP]:Envir.[Island]:Pairing[KE] |−0.05                    |0.40 |−0.13 |0.894  |
# |Dependency[RP]:Envir.[Long]:Pairing[ME]   |−0.52                    |0.40 |−1.31 |0.190  |
# |Dependency[RP]:Envir.[Island]:Pairing[ME] |−0.17                    |0.40 |−0.44 |0.661  |
# |Dependency[RP]:Envir.[Long]:Pairing[KK]   |0.08                     |0.40 |0.19  |0.847  |
# |Dependency[RP]:Envir.[Island]:Pairing[KK] |−0.63                    |0.40 |−1.56 |0.119  |
# |Dependency[RP]:Envir.[Long]:Pairing[MM]   |−0.42                    |0.39 |−1.08 |0.282  |
# |Dependency[RP]:Envir.[Island]:Pairing[MM] |−0.64                    |0.39 |−1.61 |0.107  |

# joint tests ...

joint_tests(model1)

# model term                     df1     df2 F.ratio p.value
# dependency                       1  190.14  12.497  0.0005 ***
# environment                      2  188.94  21.298  <.0001 ***
# pairing                          4  342.41  10.051  <.0001 ***
# dependency:environment           2 1251.18  11.752  <.0001 ***
# dependency:pairing               4  515.55   2.297  0.0580
# environment:pairing              8  477.80   1.190  0.3033
# dependency:environment:pairing   8 1252.59   1.411  0.1873

# post-hoc tests ...

pairwise <- model1 %>%
  emmeans(~ dependency * environment * pairing) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment pairing estimate    SE   df t.ratio p.value
# pronoun - gap short       EE        0.1667 0.208 1373   0.800  0.9997
# pronoun - gap long        EE       -0.2446 0.209 1380  -1.169  0.9842
# pronoun - gap island      EE        0.1598 0.209 1380   0.764  0.9999
# pronoun - gap short       KE        0.2273 0.199 1373   1.145  0.9871
# pronoun - gap long        KE       -0.2576 0.199 1373  -1.297  0.9606
# pronoun - gap island      KE        0.1667 0.199 1373   0.840  0.9995
# pronoun - gap short       ME        0.4852 0.190 1374   2.552  0.1501
# pronoun - gap long        ME       -0.4454 0.190 1374  -2.343  0.2526
# pronoun - gap island      ME        0.3046 0.190 1374   1.603  0.8225
# pronoun - gap short       KK        0.6970 0.199 1373   3.511  0.0069 **
# pronoun - gap long        KK        0.3636 0.199 1373   1.832  0.6465
# pronoun - gap island      KK        0.0606 0.199 1373   0.305  1.0000
# pronoun - gap short       MM        0.8630 0.189 1373   4.572  0.0001
# pronoun - gap long        MM        0.0274 0.189 1373   0.145  1.0000
# pronoun - gap island      MM        0.2192 0.189 1373   1.161  0.9852

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_iqr_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_iqr_posthoc.csv')
dt <- read_rds('tables/src_ajt_iqr_posthoc.rds')

# plot ...

plot(emmeans(model1, ~ dependency * environment | pairing), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: iqr proficiency ####
# ---------------------------------------------------------------------------- #

# remove participants with no proficiency scores ...

md <- ajt_iqr %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(proficiency = as.numeric(proficiency),
         dependency = as.factor(dependency),
         environment = as.factor(environment),
         pairing = as.factor(pairing)) %>%
  mutate(pairing = fct_drop(pairing),
         environment = fct_relevel(environment, 'short', 'long', 'island'))

# standardize proficiency scores ...

temp <- md %>%
  group_by(study, participant, proficiency) %>%
  summarise() %>%
  ungroup() %>%
  group_by(study) %>%
  mutate(proficiency = scale(proficiency))

md <- md %>%
  select(-proficiency) %>%
  left_join(temp, by = c('study', 'participant'))

summary(md)

# set contrasts ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol = 2)
contrasts(md$pairing) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol = 1)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$pairing)

# plot ...

plot <- md %>%
  mutate(group = case_when(pairing == 'KE' ~ 'KLE', 
                           pairing == 'ME' ~ 'MLE')) %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in English<br>(*n* = ', pairing_size, ')')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = proficiency, y = iqr, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Proficiency (*z*-score)'),
  scale_y_continuous(name='IQR'),
  #coord_cartesian(ylim = c(0.9, 6.1)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(-1, 3))
ggsave('plots/orc_ajt_iqr_proficiency.png', width=6.5, height=3, dpi=600)

p2 + s + coord_cartesian(ylim = c(-1, 3))
ggsave('plots/src_ajt_iqr_proficiency.png', width=6.5, height=3, dpi=600)

# ---------------------------------------------------------------------------- #
# model: iqr proficiency ####
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic() 
model2 <- lmer(iqr ~ dependency * environment * pairing * proficiency + 
                 (1 + dependency | participant),
               data = filter(md, study == 'ORC'), 
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model2)
toc()
beep()

# 0.86 sec

performance::check_convergence(model2) # okay
performance::check_singularity(model2) # okay 

# save ...

model2 %>% write_rds('models/orc_ajt_iqr_proficiency_model.rds')
model2 <- read_rds('models/orc_ajt_iqr_proficiency_model.rds')

# model summary ...

summary(model2)

# Parameter                                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                    8.899e-01  5.744e-02  1.400e+02  15.492  < 2e-16 ***
# dependency2                                   -7.293e-03  1.059e-01  1.408e+02  -0.069   0.9452    
# environment2                                   1.608e-01  9.807e-02  3.372e+02   1.640   0.1019    
# environment3                                   2.590e-01  1.175e-01  1.834e+02   2.204   0.0287 *  
# pairing2                                       1.445e-01  1.149e-01  1.400e+02   1.258   0.2104    
# proficiency                                   -2.412e-01  5.800e-02  1.400e+02  -4.160 5.53e-05 ***
# dependency2:environment2                       2.949e-01  1.929e-01  4.200e+02   1.529   0.1270    
# dependency2:environment3                       2.522e-01  1.929e-01  4.200e+02   1.308   0.1917    
# dependency2:pairing2                           8.766e-04  2.118e-01  1.408e+02   0.004   0.9967    
# environment2:pairing2                         -3.836e-02  1.961e-01  3.372e+02  -0.196   0.8451    
# environment3:pairing2                         -2.437e-02  2.350e-01  1.834e+02  -0.104   0.9175    
# dependency2:proficiency                        2.235e-01  1.069e-01  1.408e+02   2.090   0.0384 *  
# environment2:proficiency                      -1.166e-01  9.902e-02  3.372e+02  -1.177   0.2400    
# environment3:proficiency                       1.057e-01  1.186e-01  1.834e+02   0.891   0.3740    
# pairing2:proficiency                           8.897e-02  1.160e-01  1.400e+02   0.767   0.4443    
# dependency2:environment2:pairing2              4.139e-01  3.858e-01  4.200e+02   1.073   0.2839    
# dependency2:environment3:pairing2              8.831e-02  3.858e-01  4.200e+02   0.229   0.8190    
# dependency2:environment2:proficiency           1.372e-01  1.948e-01  4.200e+02   0.705   0.4814    
# dependency2:environment3:proficiency          -2.076e-01  1.948e-01  4.200e+02  -1.066   0.2871    
# dependency2:pairing2:proficiency              -3.112e-02  2.139e-01  1.408e+02  -0.146   0.8845    
# environment2:pairing2:proficiency             -7.726e-02  1.980e-01  3.372e+02  -0.390   0.6967    
# environment3:pairing2:proficiency             -2.004e-01  2.373e-01  1.834e+02  -0.845   0.3994    
# dependency2:environment2:pairing2:proficiency -8.573e-02  3.895e-01  4.200e+02  -0.220   0.8259    
# dependency2:environment3:pairing2:proficiency  4.486e-01  3.895e-01  4.200e+02   1.152   0.2501 

# table ...

dt <- summary(model2)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_iqr_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_iqr_proficiency_model.csv')
dt <- read_rds('tables/orc_ajt_iqr_proficiency_model.rds')

# joint tests ...

joint_tests(model2)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 140   0.005  0.9451
# environment                                  2 560   3.131  0.0444 *
# pairing                                      1 140   1.583  0.2104
# proficiency                                  1 140  17.303  0.0001 ***
# dependency:environment                       2 560   1.163  0.3133
# dependency:pairing                           1 140   0.000  0.9967
# dependency:proficiency                       1 140   4.378  0.0382 *
# environment:pairing                          2 560   0.017  0.9829
# environment:proficiency                      2 560   2.221  0.1095
# pairing:proficiency                          1 140   0.588  0.4443
# dependency:environment:pairing               2 560   0.544  0.5807
# dependency:environment:proficiency           2 560   1.354  0.2592
# dependency:pairing:proficiency               1 140   0.021  0.8844
# environment:pairing:proficiency              2 560   0.459  0.6322
# dependency:environment:pairing:proficiency   2 560   0.925  0.3973

pairwise <- model2 %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment pairing estimate    SE  df t.ratio p.value
# pronoun - gap short       KE         0.323 0.243 521   1.330  0.6985
# pronoun - gap long        KE         0.503 0.243 521   2.071  0.2090
# pronoun - gap island      KE        -0.109 0.243 521  -0.448  0.9982
# pronoun - gap short       ME         0.171 0.214 521   0.797  0.9624
# pronoun - gap long        ME         0.265 0.214 521   1.237  0.7627
# pronoun - gap island      ME         0.188 0.214 521   0.875  0.9417

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/orc_ajt_iqr_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ajt_iqr_proficiency_posthoc.csv')
dt <- read_rds('tables/orc_ajt_iqr_proficiency_posthoc.rds')

# plot ... 

emmip(model2, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

tic() 
model3 <- lmer(iqr ~ dependency * environment * pairing * proficiency + 
                 (1 | participant),
               data = filter(md, study == 'SRC'), 
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model3)
toc()
beep()

# 1 sec

performance::check_convergence(model3) # okay
performance::check_singularity(model3) # okay

# save ...

model3 %>% write_rds('models/src_ajt_iqr_proficiency_model.rds')
model3 <- read_rds('models/src_ajt_iqr_proficiency_model.rds')

# model summary ...

summary(model3)

# Parameter                                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     0.879784   0.062868 133.000000  13.994  < 2e-16 ***
# dependency2                                     0.095276   0.084236 665.000000   1.131 0.258436    
# environment2                                    0.433475   0.103167 665.000000   4.202 3.01e-05 ***
# environment3                                    0.401612   0.103167 665.000000   3.893 0.000109 ***
# pairing2                                        0.053422   0.125737 133.000000   0.425 0.671613    
# proficiency                                    -0.075617   0.063229 133.000000  -1.196 0.233861    
# dependency2:environment2                       -0.874306   0.206335 665.000000  -4.237 2.58e-05 ***
# dependency2:environment3                       -0.113768   0.206335 665.000000  -0.551 0.581561    
# dependency2:pairing2                            0.145903   0.168471 665.000000   0.866 0.386781    
# environment2:pairing2                          -0.365805   0.206335 665.000000  -1.773 0.076707 .  
# environment3:pairing2                          -0.539035   0.206335 665.000000  -2.612 0.009193 ** 
# dependency2:proficiency                         0.148148   0.084720 665.000000   1.749 0.080808 .  
# environment2:proficiency                       -0.112561   0.103760 665.000000  -1.085 0.278394    
# environment3:proficiency                       -0.035667   0.103760 665.000000  -0.344 0.731147    
# pairing2:proficiency                            0.233461   0.126459 133.000000   1.846 0.067095 .  
# dependency2:environment2:pairing2              -0.521233   0.412669 665.000000  -1.263 0.207004    
# dependency2:environment3:pairing2              -0.125527   0.412669 665.000000  -0.304 0.761084    
# dependency2:environment2:proficiency           -0.080850   0.207520 665.000000  -0.390 0.696958    
# dependency2:environment3:proficiency           -0.045465   0.207520 665.000000  -0.219 0.826650    
# dependency2:pairing2:proficiency                0.058486   0.169439 665.000000   0.345 0.730072    
# environment2:pairing2:proficiency               0.192506   0.207520 665.000000   0.928 0.353925    
# environment3:pairing2:proficiency              -0.025452   0.207520 665.000000  -0.123 0.902422    
# dependency2:environment2:pairing2:proficiency  -1.486362   0.415040 665.000000  -3.581 0.000367 ***
# dependency2:environment3:pairing2:proficiency   0.007786   0.415040 665.000000   0.019 0.985039  

# table ...

dt <- summary(model3)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  mutate_if(is.character, str_replace_all, 
            c('pairing2' = 'pairing3')) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_iqr_proficiency_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_iqr_proficiency_model.csv')
dt <- read_rds('tables/src_ajt_iqr_proficiency_model.rds')

# joint tests ...

joint_tests(model3)

# model term                                 df1 df2 F.ratio p.value
# dependency                                   1 665   1.279  0.2584
# environment                                  2 665  10.968  <.0001 ***
# pairing                                      1 133   0.181  0.6716
# proficiency                                  1 133   1.430  0.2339
# dependency:environment                       2 665  10.615  <.0001 ***
# dependency:pairing                           1 665   0.750  0.3868
# dependency:proficiency                       1 665   3.058  0.0808 .
# environment:pairing                          2 665   3.558  0.0291 *
# environment:proficiency                      2 665   0.615  0.5411
# pairing:proficiency                          1 133   3.408  0.0671 .
# dependency:environment:pairing               2 665   0.869  0.4198
# dependency:environment:proficiency           2 665   0.076  0.9266
# dependency:pairing:proficiency               1 665   0.119  0.7301
# environment:pairing:proficiency              2 665   0.660  0.5174
# dependency:environment:pairing:proficiency   2 665   8.595  0.0002 ***

pairwise <- model3 %>%
  emtrends(~ dependency * environment * pairing, 
           var = 'proficiency') %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# post-hoc tests ...

# contrast      environment pairing estimate    SE  df t.ratio p.value
# pronoun - gap short       KE       -0.0854 0.215 665  -0.397  0.9991
# pronoun - gap long        KE        0.5769 0.215 665   2.679  0.0445 *
# pronoun - gap island      KE       -0.1348 0.215 665  -0.626  0.9894
# pronoun - gap short       ME        0.4659 0.199 665   2.337  0.1126
# pronoun - gap long        ME       -0.3581 0.199 665  -1.796  0.3645
# pronoun - gap island      ME        0.4244 0.199 665   2.128  0.1855

# table ...

dt <- pairwise %>%
  select(pairing, environment, contrast, everything()) %>%
  select(-df) %>%
  arrange(pairing, environment) %>%
  my_table()

dt %>% write_rds('tables/src_ajt_iqr_proficiency_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ajt_iqr_proficiency_posthoc.csv')
dt <- read_rds('tables/src_ajt_iqr_proficiency_posthoc.rds')

# plot ... 

emmip(model3, dependency ~ proficiency | pairing + environment, cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: saturation ####
# ---------------------------------------------------------------------------- #

# remove participants with no proficiency scores ...

ajt_index <- ajt_clean %>%
  filter(dependency %in% c('gap', 'pronoun'),
         environment %in% c('short', 'long', 'island')) %>%
  arrange(study, pairing, participant, condition, trial_index) %>%
  group_by(study, pairing, participant, condition) %>%
  mutate(index = row_number()) %>%
  ungroup() %>%
  mutate(dependency = as.factor(dependency),
         environment = as.factor(environment),
         pairing = as.factor(pairing)) %>%
  mutate(pairing = fct_drop(pairing),
         environment = fct_relevel(environment, 'short', 'long', 'island'))

# plot 1 ...

plot <- ajt_index %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(*n* = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, dependency, environment, index) %>%
  summarise(mean = mean(response, na.rm=T),
            lcl = my.lcl(response),
            ucl = my.ucl(response)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 
                                   'short' = 'Short', 
                                   'long' = 'Long',
                                   'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

design <- '
ABD
#CE
'

s <- list(
  aes(x=index, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency, group = dependency), position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=.3, linetype=1, 
                position = position_dodge(width = .4),
                show.legend = FALSE),
  geom_point(aes(col = dependency), size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(1, 6)),
  scale_x_continuous(name="Exposure"),
  scale_y_continuous(name="Rating",
                     breaks = c(1, 3.5, 6)),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c("Gap", "RP")),
  scale_shape_manual(name="Dependency", 
                     values=c(16, 15), labels=c("Gap", "RP")),
  scale_linetype_manual(name="Dependency", 
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(
    text = element_text(size = 12),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = 'bottom',
    legend.margin = margin(-10, 0, 0, 0),
    strip.text.x = element_markdown(),
    strip.text.y = element_markdown()),
  facet_grid2(pairing_label~environment, 
               axes = 'all', 
               remove_labels = FALSE) 
)

p1 + s
ggsave("plots/orc_ajt_exposure1.png", width=6.5, height=7, dpi=600)

p2 + s
ggsave("plots/src_ajt_exposure1.png", width=6.5, height=7, dpi=600)

# plot 2 ...

plot <- ajt_index %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '<br>(*n* = ', pairing_size, ')')) %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = index, y = response, col = dependency, fill = dependency),
  geom_smooth(aes(linetype = dependency), linewidth = .5, alpha = .75, method=lm),
  theme_classic(),
  scale_x_continuous(name='Exposure'),
  scale_y_continuous(name='Rating',
                     limits = c(1, 6),
                     breaks = c(1, 3.5, 6)),
  scale_color_manual(name = 'Dependency',
                     values=c('#648fff', '#ffb000'), 
                     labels=c("Gap", "RP")),
  scale_fill_manual(name = 'Dependency',
                    values=c('#c1d2ff', '#ffdf99'), 
                    labels=c("Gap", "RP")),
  scale_linetype_manual(name = 'Dependency',
                        values=c('41', 'solid'), 
                        labels=c("Gap", "RP")),
  theme(text = element_text(size = 12),
        strip.text.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "right"),
  facet_grid2(pairing_label~environment, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ajt_exposure2.png', width=6.5, height=7.5, dpi=600)

p2 + s
ggsave('plots/src_ajt_exposure2.png', width=6.5, height=7.5, dpi=600)

# ---------------------------------------------------------------------------- #
# plot: RP acceptors vs. non-acceptors (bar plot) ####
# ---------------------------------------------------------------------------- #

# plot by group ...

ajt_ind2 <- ajt_ind %>%
  group_by(study, group, language, pairing, participant, proficiency) %>%
  summarise(max = max(acceptance)) %>%
  ungroup() %>%
  mutate(subgroup = case_when(max >= .8 ~ "RP Acceptor",
                              max <= .2 ~ "RP Rejector",
                              TRUE ~ "Inconsistent"))

plot <- ajt_ind2 %>%
  mutate(subgroup = case_when(subgroup == "RP Acceptor" ~ "RP Acceptor",
                              TRUE ~ "RP Non-Acceptor")) %>%
  filter(language == "English") %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study, group, group_label, subgroup) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(study, group) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(prop = count / total) %>%
  mutate(subgroup = fct_relevel(subgroup, 
                                "RP Acceptor", 
                                "RP Non-Acceptor"
  ))

p <- ggplot(data=plot)
p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

# figure without plus signs and minus signs ...

s <- list(
  aes(x=group, y=(count/total)*100, group=subgroup, fill=subgroup, col=subgroup, 
      label=paste0(round(count/total * 100), '% (*n* = ', count, ')')),
  geom_bar(stat = "identity", width = .9, 
           col = 'black', 
           #linewidth = .25
  ),
  geom_richtext(col = 'black', fill = NA, label.color = NA, size = 3,
                position = position_stack(vjust = 0.5)),
  theme_classic(),
  scale_x_discrete(name='Group'),
  scale_y_continuous(name='Participants (%)', 
  ),
  scale_fill_manual(name="Type:", 
                    breaks = c('RP Acceptor', 'RP Non-Acceptor'),
                    #labels = c('RP Acceptor', 'RP Non-Acceptor'),
                    values=c('#CDFAD5', '#ffbfbf')),
  theme(text = element_text(size = 12),
        axis.text.x = element_markdown(),
        strip.text = element_markdown(),
        #legend.text = element_text(),
        legend.box.margin = margin(-7, 0, 0, 0),
        legend.position = "bottom", 
        legend.margin=margin(1, 1, 1, 1)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x'),
  facet_wrap(~study, scales = "free_x")
) 

p + s
ggsave('plots/ajt_subgroup_acceptor.png', width=6.5, height=2.5, dpi=600)

# figure with plus signs and minus signs ...

s <- list(
  aes(x=group, y=(count/total)*100, group=subgroup, fill=subgroup, col=subgroup, 
      label=paste0(round(count/total * 100), '% (*n* = ', count, ')')),
  geom_bar(stat = "identity", width = .9, 
           col = 'black', 
           #linewidth = .25
  ),
  geom_richtext(col = 'black', fill = NA, label.color = NA, size = 3,
                position = position_stack(vjust = 0.5)),
  theme_classic(),
  scale_x_discrete(name='Group'),
  scale_y_continuous(name='Participants (%)', 
  ),
  scale_fill_manual(name="Type:", 
                    breaks = c('RP Acceptor', 'RP Non-Acceptor'),
                    labels = c('+ RP Acceptor', '\u2212 RP Acceptor'),
                    #labels = c('RP Acceptor', 'RP Non-Acceptor'),
                    values=c('#CDFAD5', '#ffbfbf')),
  theme(text = element_text(size = 12),
        axis.text.x = element_markdown(),
        strip.text = element_markdown(),
        #legend.text = element_text(),
        legend.box.margin = margin(-7, 0, 0, 0),
        legend.position = "bottom", 
        legend.margin=margin(1, 1, 1, 1)),
  facet_manual(~pairing_label, design = 'ABDCE', scales = 'free_x'),
  facet_wrap(~study, scales = "free_x")
) 

p + s
ggsave('plots/ajt_subgroup_acceptor.png', width=6.5, height=2.5, dpi=600)

p1 + s
ggsave('plots/orc_ajt_subgroup_acceptor.png', width=6.5, height=2, dpi=600)

p2 + s
ggsave('plots/src_ajt_subgroup_acceptor.png', width=6.5, height=2, dpi=600)

acceptor <- ajt_ind2 %>%
  filter(subgroup == "RP Acceptor", language == "English") %>%
  select(-max) 

check <- acceptor %>%
  group_by(study, pairing) %>%
  summarise(n = length(unique(participant))) %>%
  ungroup()

acceptor %>% write_csv("data/acceptor.csv")

# ============================================================================ #
# ** REANALYSIS OF EPT AND SPRT DATA FOR RP NON-ACCEPTORS ** ####
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# plot: AJT (RP acceptance) vs. EPT (token) ####
# ---------------------------------------------------------------------------- #

ajt_ind2 <- ajt_ind %>%
  group_by(study, group, language, pairing, participant, proficiency) %>%
  summarise(max = max(acceptance)) %>%
  ungroup() %>%
  mutate(subgroup = case_when(max >= .8 ~ "RP Acceptor",
                          max <= .2 ~ "RP Rejector",
                          TRUE ~ "Inconsistent")) %>%
  filter(language == 'English')

check <- ajt_ind2 %>%
  group_by(study, pairing) %>%
  summarise(n = n()) %>%
  ungroup()

ep2 <- ep %>%
  left_join(ajt_ind2, by = c('study', 'group', 'participant')) %>%
  filter(subgroup %in% c('RP Acceptor', 'RP Rejector', 'Inconsistent')) %>%
  mutate(subgroup = fct_relevel(subgroup, 'RP Acceptor', 'RP Rejector', 'Inconsistent')) %>%
  mutate(subgroup2 = case_when(subgroup %in% c('RP Rejector', 'Inconsistent') ~ 'RP Non-Acceptor',
                               subgroup == 'RP Acceptor' ~ 'RP Acceptor')) %>%
  mutate(subgroup2 = fct_relevel(subgroup2, 'RP Acceptor', 'RP Non-Acceptor'))

# summarise for plotting ...

plot <- ep2 %>%
  filter(is.na(condition) == FALSE) %>%
  mutate(type = as.character(type)) %>%
  group_by(study, group, subgroup) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '<br>(*n* = ', group_size, ')')) %>%
  group_by(study, group, subgroup, group_size, group_label, participant, condition, type) %>%
  count(type) %>%
  ungroup() %>%
  complete(type, nesting(study, group, subgroup, group_size, group_label, participant, condition), fill = list(n = 0)) %>%
  select(study, group, subgroup, group_size, group_label, participant, condition, type, n) %>%
  arrange(study, group, subgroup, group_size, group_label, participant, condition) %>%
  group_by(study, group, subgroup, group_size, group_label, condition, type) %>%
  summarise(
    mean = mean(n, na.rm = TRUE),
    lcl = my.lcl(n, conf.level = .95),
    ucl = my.ucl(n, conf.level = .95)
  ) %>%
  mutate(lcl = case_when(lcl == 'NaN' ~ 0, TRUE ~ lcl),
         ucl = case_when(ucl == 'NaN' ~ 0, TRUE ~ ucl),
         prop = mean * 20
  ) %>%
  mutate(lcl = case_when(is.na(lcl) == TRUE ~ mean, TRUE ~ lcl),
         ucl = case_when(is.na(ucl) == TRUE ~ mean, TRUE ~ ucl)) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = c('gap', 'other', 'nontarget', 'resumption'))) %>%
  filter(!type %in% c('nontarget'))

plot2 <- ep2 %>%
  mutate(subgroup = subgroup2) %>%
  filter(is.na(condition) == FALSE) %>%
  mutate(type = as.character(type)) %>%
  group_by(study, group, subgroup) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, ' (*n* = ', group_size, ')')) %>%
  mutate(group_label2 = paste0(group, ' (*n* = ', group_size, ')')) %>%
  group_by(study, group, subgroup, group_size, group_label, group_label2, participant, condition, type) %>%
  count(type) %>%
  ungroup() %>%
  complete(type, nesting(study, group, subgroup, group_size, group_label, group_label2, participant, condition), fill = list(n = 0)) %>%
  select(study, group, subgroup, group_size, group_label, group_label2, participant, condition, type, n) %>%
  arrange(study, group, subgroup, group_size, group_label, group_label2, participant, condition) %>%
  group_by(study, group, subgroup, group_size, group_label, group_label2, condition, type) %>%
  summarise(
    mean = mean(n, na.rm = TRUE),
    lcl = my.lcl(n, conf.level = .95),
    ucl = my.ucl(n, conf.level = .95)
  ) %>%
  mutate(lcl = case_when(lcl == 'NaN' ~ 0, TRUE ~ lcl),
         ucl = case_when(ucl == 'NaN' ~ 0, TRUE ~ ucl),
         prop = mean * 20
  ) %>%
  mutate(lcl = case_when(is.na(lcl) == TRUE ~ mean, TRUE ~ lcl),
         ucl = case_when(is.na(ucl) == TRUE ~ mean, TRUE ~ ucl)) %>%
  ungroup() %>%
  mutate(type = factor(type, levels = c('gap', 'other', 'nontarget', 'resumption'))) %>%
  filter(!type %in% c('nontarget')) %>%
  mutate(subgroup2 = case_when(subgroup == 'RP Acceptor' ~ '+ RP Acceptor', 
                               TRUE ~ '\u2212 RP Acceptor'))

# define data for plots ...

p3a <- ggplot(data=filter(plot2, study == 'ORC', subgroup == 'RP Acceptor'))
p3b <- ggplot(data=filter(plot2, study == 'ORC', subgroup == 'RP Non-Acceptor'))
p4a <- ggplot(data=filter(plot2, study == 'SRC', subgroup == 'RP Acceptor'))
p4b <- ggplot(data=filter(plot2, study == 'SRC', subgroup == 'RP Non-Acceptor'))

# just non-acceptors ...

# styling ...

s <- list(
  aes(x=condition, y=mean, group=type, col=type, shape=type, label=paste0(round(prop), '%')),
  geom_hline(yintercept = 2.5),
  geom_line(aes(linetype = type), position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1,
             shape = 21,
             fill = 'white', 
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(0, 5)),
  scale_y_continuous(name='Tokens', breaks = c(0, 2.5, 5)),
  scale_x_discrete(name='Environment', 
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000', 'coral'),
                      breaks = c('gap', 'resumption', 'other'),
                      labels = c('Gapping', 'Resumption', 'Other')),
  scale_linetype_manual(name='Dependency', 
                        values=c('41', 'solid', '7222'),
                        breaks = c('gap', 'resumption', 'other'),
                        labels = c('Gapping', 'Resumption', 'Other')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'right',
        #legend.position = 'bottom',
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

p3b + s
ggsave('plots/orc_ept_token_nonacceptor.png', width=6.5, height=1.75, dpi=600)

p4b + s
ggsave('plots/src_ept_token_nonacceptor.png', width=6.5, height=1.75, dpi=600)

# for defense (resumption only) ...

p5 <- ggplot(data=filter(plot2, study == 'ORC', subgroup == 'RP Non-Acceptor', type == 'resumption'))

s2 <- list(
  aes(x=condition, y=mean*20, group = type),
  geom_bar(stat = "identity", col = 'black', linewidth = .5, fill = '#ffb000', width = .7, alpha=.7),
  theme_classic(),
  scale_y_continuous(name='Responses (%)', 
                     limits = c(0, 40),
                     #breaks = c(0, 25, 50)
                     ),
  scale_x_discrete(name='Environment', 
                   labels = c('Short', 'Long', 'Island')
  ),
  theme(text = element_text(size = 12),
        plot.margin = margin(2, 2, 0, 2),
        #legend.position = 'right', # dissertation
        legend.position = 'none', # defense
        # legend.margin = margin(t = -.4, unit = 'cm'),
        #legend.key.width = unit(.8, 'cm'),
        plot.title = element_text(size = 12, hjust = .5),
        strip.text = element_markdown()),
  facet_manual(~group_label,
               design = 'ABC',
               axes = 'all',
               remove_labels = FALSE)
)

p5 + s2
ggsave('plots/defense/ept_token_nonacceptor.png', width=5, height=2.5, dpi=600)

plot3a = p3a + s + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
plot3b = p3b + s + theme(axis.title.y = element_blank())
grid.arrange(patchworkGrob(plot3a / plot3b +
                             plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                             plot_layout(guides = 'collect')), left = 'Tokens')
g <- arrangeGrob(patchworkGrob(plot3a / plot3b +
                                 plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                                 plot_layout(guides = 'collect')), left = 'Tokens')
ggsave('plots/orc_ept_token_subgroup2.png', g, width=6.5, height=3.25, dpi=600)

plot4a = p4a + s + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
plot4b = p4b + s + theme(axis.title.y = element_blank())
grid.arrange(patchworkGrob(plot4a / plot4b + 
                             plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                             plot_layout(guides = 'collect')), left = 'Tokens')
g <- arrangeGrob(patchworkGrob(plot4a / plot4b + 
                                 plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                                 plot_layout(guides = 'collect')), left = 'Tokens')
ggsave('plots/src_ept_token_subgroup2.png', g, width=6.5, height=3.25, dpi=600)

# ---------------------------------------------------------------------------- #
# model: AJT (RP acceptance) vs. EPT (token) ####
# ---------------------------------------------------------------------------- #

# prep data for orc modeling ...

temp <- ajt_ind2 %>%
  filter(subgroup %in% c('RP Rejector', 'Inconsistent'))

check <- temp %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

temp2 <- ep %>%
  filter(participant %in% temp$participant) %>%
  mutate(resumption = case_when(type == 'resumption' ~ TRUE,
                                TRUE ~ FALSE)) %>%
  filter(group %in% c('ENS', 'KLE', 'MLE')) %>%
  mutate(group = fct_drop(group)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  select(study, group, participant, item, environment, resumption, type) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'environment', 'type'), factor) %>%
  filter(is.na(environment) == FALSE) %>%
  mutate(environment = fct_drop(environment))

# add dummy participants with only resumption responses to the SRC study's ENS and KLE groups ...

dummy <- read_csv('data/ept_dummy.csv') %>% 
  filter(group %in% c('ENS', 'KLE'))

md <- bind_rows(temp2, dummy) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'environment', 'type'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(study = fct_drop(study)) %>%
  mutate(environment = fct_drop(environment))

summary(md)

check <- md %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model4 <- glmer(resumption ~ environment * group +
                  (1 | participant) + 
                  (1 + group | item), 
                data = filter(md, study == 'ORC'), family = binomial, 
                control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
summary(model4)
toc()
beep()

# ... 11.45 sec elapsed
# ... correlations of 1 ??

check_convergence(model4) # okay
check_singularity(model4) # warning (random effects correlations look okay, though)
MuMIn::r.squaredGLMM(model4) # R2m = ???, R2c = ???

# save ...

model4 %>% write_rds('models/orc_ept_token_subgroup_model.rds')
model4 <- read_rds('models/orc_ept_token_subgroup_model.rds') 

# model summary ...

summary(model4)

# Parameter           Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -5.4837     0.5027 -10.909  < 2e-16 ***
# environment2          1.8623     0.5752   3.238 0.001206 ** 
# environment3          3.9060     0.5928   6.589 4.43e-11 ***
# group2               -1.2409     0.9137  -1.358 0.174421    
# group3               -0.2556     0.8259  -0.309 0.756998    
# environment2:group2   1.1773     1.5339   0.768 0.442776    
# environment3:group2  -0.6863     1.5479  -0.443 0.657486    
# environment2:group3  -0.1051     1.0733  -0.098 0.921963    
# environment3:group3  -3.9421     1.0597  -3.720 0.000199 ***

# table ...

dt <- summary(model4)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/orc_ept_token_subgroup_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ept_token_subgroup_model.csv')
dt <- read_rds('tables/orc_ept_token_subgroup_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                      |Estimate~($\hat{\beta}$) |$SE$ |$z$    |$p$    |
# |:------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                    |−5.48                    |0.50 |−10.91 |< .001 |
# |Environment[Long]              |1.86                     |0.58 |3.24   |0.001  |
# |Environment[Island]            |3.91                     |0.59 |6.59   |< .001 |
# |Group[KLE]                     |−1.24                    |0.91 |−1.36  |0.174  |
# |Group[MLE]                     |−0.26                    |0.83 |−0.31  |0.757  |
# |Environment[Long]:Group[KLE]   |1.18                     |1.53 |0.77   |0.443  |
# |Environment[Island]:Group[KLE] |−0.69                    |1.55 |−0.44  |0.657  |
# |Environment[Long]:Group[MLE]   |−0.11                    |1.07 |−0.10  |0.922  |
# |Environment[Island]:Group[MLE] |−3.94                    |1.06 |−3.72  |< .001 |

# joint tests ...

joint_tests(model4)

# model term        df1 df2 F.ratio p.value
# environment         2 Inf  28.483  <.0001
# group               2 Inf   0.925  0.3966
# environment:group   4 Inf   7.904  <.0001

# post-hoc tests ...

pairwise <- model4 %>% 
  emmeans(~ environment * group) %>%
  contrast('revpairwise', simple = 'environment', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast       group estimate    SE  df z.ratio p.value
# long - short   ENS      1.505 0.737 Inf   2.042  0.2675
# island - short ENS      5.449 0.738 Inf   7.384  <.0001 ***
# island - long  ENS      3.944 0.468 Inf   8.436  <.0001 ***
# long - short   KLE      2.682 1.344 Inf   1.995  0.2936
# island - short KLE      4.762 1.385 Inf   3.439  0.0048 **
# island - long  KLE      2.080 0.746 Inf   2.787  0.0415 *
# long - short   MLE      1.400 0.784 Inf   1.786  0.4257
# island - short MLE      1.507 0.774 Inf   1.947  0.3219
# island - long  MLE      0.107 0.587 Inf   0.182  1.0000

# table ...

dt <- pairwise %>%
  select(group, contrast, estimate, SE, z.ratio, p.value) %>%
  my_table()

dt %>% write_rds('tables/orc_ept_token_subgroup_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_ept_token_subgroup_posthoc.csv')
dt <- read_rds('tables/orc_ept_token_subgroup_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Contrast       |Estimate~($\hat{\beta}$) |$SE$ |$z$  |$p$    |
# |:-----|:--------------|:------------------------|:----|:----|:------|
# |ENS   |Long − Short   |1.50                     |0.74 |2.04 |0.267  |
# |ENS   |Island − Short |5.45                     |0.74 |7.38 |< .001 |
# |ENS   |Island − Long  |3.94                     |0.47 |8.44 |< .001 |
# |KLE   |Long − Short   |2.68                     |1.34 |2.00 |0.294  |
# |KLE   |Island − Short |4.76                     |1.38 |3.44 |0.005  |
# |KLE   |Island − Long  |2.08                     |0.75 |2.79 |0.042  |
# |MLE   |Long − Short   |1.40                     |0.78 |1.79 |0.426  |
# |MLE   |Island − Short |1.51                     |0.77 |1.95 |0.322  |
# |MLE   |Island − Long  |0.11                     |0.59 |0.18 |> .999 |

# plot ...

plot(emmeans(model4, ~ environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model2 <- glmer(resumption ~ environment * group +
                  (1 | participant) + 
                  (1 | item), 
                data = filter(md, study == 'SRC'), family = binomial, 
                control = glmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 1e6)))
summary(model2)
toc()
beep()

# 9 sec
# no warnings

check_convergence(model2) # okay
check_singularity(model2) # okay

# set as final model ...

model <- model2

# save ...

model %>% write_rds('models/src_ept_token_subgroup_model.rds')
model <- read_rds('models/src_ept_token_subgroup_model.rds') 

# model summary ...

summary(model)

# table ...

dt <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_ept_token_subgroup_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ept_token_subgroup_model.csv')
dt <- read_rds('tables/src_ept_token_subgroup_model.rds')

dt %>% rename('Est.' = 2) %>% kbl(format = 'markdown')

# |Parameter                      |Est.  |$SE$ |$z$   |$p$    |
# |:------------------------------|:-----|:----|:-----|:------|
# |(Intercept)                    |−8.23 |1.57 |−5.25 |< .001 | ***
# |Environment[Long]              |6.16  |1.30 |4.73  |< .001 | ***
# |Environment[Island]            |11.37 |1.95 |5.83  |< .001 | ***
# |Group[KLE]                     |−4.82 |2.91 |−1.66 |0.098  | .
# |Group[MLE]                     |−4.01 |2.79 |−1.44 |0.150  |
# |Environment[Long]:Group[KLE]   |6.22  |2.94 |2.12  |0.034  | *
# |Environment[Island]:Group[KLE] |−6.26 |3.62 |−1.73 |0.084  | .
# |Environment[Long]:Group[MLE]   |5.11  |2.64 |1.93  |0.053  | .
# |Environment[Island]:Group[MLE] |−9.97 |3.93 |−2.54 |0.011  | *

# joint tests ...

joint_tests(model)

# model term        df1 df2 F.ratio p.value
# environment         2 Inf  17.098  <.0001 ***
# group               2 Inf   1.444  0.2361
# environment:group   4 Inf   6.341  <.0001 ***

# post-hoc tests ...

pairwise <- model %>% 
  emmeans(~ environment * group) %>%
  contrast('revpairwise', simple = 'environment', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# table ...

dt <- pairwise %>%
  select(group, contrast, estimate, SE, z.ratio, p.value) %>%
  my_table()

dt %>% write_rds('tables/src_ept_token_subgroup_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_ept_token_subgroup_posthoc.csv')
dt <- read_rds('tables/src_ept_token_subgroup_posthoc.rds')

dt %>% rename('Est.' = 3) %>% kbl(format = 'markdown')

# |Group |Contrast       |Est.  |$SE$ |$z$   |$p$    |
# |:-----|:--------------|:-----|:----|:-----|:------|
# |ENS   |Long − Short   |2.38  |2.11 |1.13  |0.846  |
# |ENS   |Island − Short |16.78 |3.87 |4.34  |< .001 | ***
# |ENS   |Island − Long  |14.40 |3.00 |4.80  |< .001 | ***
# |KLE   |Long − Short   |8.60  |2.29 |3.75  |0.001  | **
# |KLE   |Island − Short |10.52 |2.40 |4.39  |< .001 | ***
# |KLE   |Island − Long  |1.92  |0.59 |3.26  |0.008  | **
# |MLE   |Long − Short   |7.49  |1.73 |4.32  |< .001 | ***
# |MLE   |Island − Short |6.81  |1.66 |4.09  |< .001 | ***
# |MLE   |Island − Long  |−0.68 |0.84 |−0.80 |0.965  |

# plot ...

plot(emmeans(model, ~ environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: AJT (RP acceptance) vs. SPR (Log RT) ####
# ---------------------------------------------------------------------------- #

ajt_ind2 <- ajt_ind %>%
  group_by(study, group, language, pairing, participant, proficiency) %>%
  summarise(max = max(acceptance)) %>%
  ungroup() %>%
  mutate(subgroup = case_when(max >= .8 ~ "RP Acceptor",
                              max <= .2 ~ "RP Rejector",
                              TRUE ~ "Inconsistent")) %>%
  filter(language == 'English')

check <- ajt_ind2 %>%
  group_by(study, pairing) %>%
  summarise(n = n()) %>%
  ungroup()

sp2 <- spr_rt %>% 
  select(study, group, participant, item, dependency, environment, region2, rt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(logrt = log(as.numeric(rt))) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  left_join(ajt_ind2, by = c('study', 'group', 'participant')) %>%
  filter(subgroup %in% c('RP Acceptor', 'RP Rejector', 'Inconsistent')) %>%
  mutate(subgroup = fct_relevel(subgroup, 'RP Acceptor', 'RP Rejector', 'Inconsistent')) %>%
  mutate(subgroup2 = case_when(subgroup %in% c('RP Rejector', 'Inconsistent') ~ 'RP Non-Acceptor',
                               subgroup == 'RP Acceptor' ~ 'RP Acceptor')) %>%
  mutate(subgroup2 = fct_relevel(subgroup2, 'RP Acceptor', 'RP Non-Acceptor'))

plot <- sp2 %>%
  group_by(study, group, subgroup) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, ' (*n* = ', group_size, ')')) %>%
  group_by(study, group, subgroup, group_size, group_label, dependency, environment, region2) %>%
  summarise(
    mean = mean(logrt, na.rm = TRUE),
    lcl = my.lcl(logrt),
    ucl = my.ucl(logrt))

plot2 <- sp2 %>%
  mutate(subgroup = subgroup2) %>%
  group_by(study, group, subgroup) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, ' (*n* = ', group_size, ')')) %>%
  group_by(study, group, subgroup, group_size, group_label, dependency, environment, region2) %>%
  summarise(
    mean = mean(logrt, na.rm = TRUE),
    lcl = my.lcl(logrt),
    ucl = my.ucl(logrt)) %>%
  ungroup() %>%
  mutate(subgroup2 = case_when(subgroup == "RP Acceptor" ~ "+ RP Acceptor",
                               TRUE ~ "\u2212 RP Acceptor"))

# define data for plots ...

p3a <- ggplot(data=filter(plot2, study == 'ORC', region2 == 1, subgroup == 'RP Acceptor'))
p3b <- ggplot(data=filter(plot2, study == 'ORC', region2 == 1, subgroup == 'RP Non-Acceptor'))
p4a <- ggplot(data=filter(plot2, study == 'SRC', region2 == 2, subgroup == 'RP Acceptor'))
p4b <- ggplot(data=filter(plot2, study == 'SRC', region2 == 2, subgroup == 'RP Non-Acceptor'))

# styling for just non-acceptors ...

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_hline(yintercept = 0),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(5.45, 6.55)),
  scale_y_continuous(name='Log RT', 
                     breaks = c(5.5, 6.0, 6.5)
  ),
  scale_x_discrete(name='Environment', 
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000'),
                      breaks = c('gap', 'pronoun'),
                      labels = c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency', 
                        values=c('41', 'solid'),
                        breaks = c('gap', 'pronoun'),
                        labels = c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'right',
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

p3b + s 
ggsave('plots/orc_spr_rt_nonacceptor.png', width=6.5, height=1.75, dpi=600)

p4b + s
ggsave('plots/src_spr_rt_nonacceptor.png', width=6.5, height=1.75, dpi=600)

# styling ...

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_hline(yintercept = 0),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1, shape = 21, fill = 'white',
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(5.45, 6.55)),
  scale_y_continuous(name='Log RT', 
                     breaks = c(5.5, 6.0, 6.5)
                     ),
  scale_x_discrete(name='Environment', 
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency', 
                      values=c('#648fff', '#ffb000'),
                      breaks = c('gap', 'pronoun'),
                      labels = c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency', 
                     values=c('41', 'solid'),
                     breaks = c('gap', 'pronoun'),
                     labels = c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'right',
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(subgroup2~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

plot3a = p3a + s + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
plot3b = p3b + s + theme(axis.title.y = element_blank())
grid.arrange(patchworkGrob(plot3a / plot3b + 
                             plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                             plot_layout(guides = 'collect')), left = 'Log RT')
g <- arrangeGrob(patchworkGrob(plot3a / plot3b + 
                                 plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                                 plot_layout(guides = 'collect')), left = 'Log RT')
ggsave('plots/orc_spr_rt_subgroup2.png', g, width=6.5, height=3.25, dpi=600)

plot4a = p4a + s + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
plot4b = p4b + s + theme(axis.title.y = element_blank())
grid.arrange(patchworkGrob(plot4a / plot4b + 
                             plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                             plot_layout(guides = 'collect')), left = 'Log RT')
g <- arrangeGrob(patchworkGrob(plot4a / plot4b + 
                                 plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                                 plot_layout(guides = 'collect')), left = 'Log RT')
ggsave('plots/src_spr_rt_subgroup2.png', g, width=6.5, height=3.25, dpi=600)

# nested plots ...

p5 <- ggplot(data=filter(plot2, study == 'ORC', region2 == '1'))
p6 <- ggplot(data=filter(plot2, study == 'SRC', region2 == '2'))

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1,
             shape = 21,
             fill = 'white',
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(5.45, 6.55)),
  scale_y_continuous(name='Log RT', breaks = c(5.5, 6, 6.5)),
  scale_x_discrete(name='Environment',
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency:',
                      values=c('#648fff', '#ffb000'),
                      breaks = c('gap', 'pronoun'),
                      labels = c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:',
                        values=c('41', 'solid'),
                        breaks = c('gap', 'pronoun'),
                        labels = c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom',
        strip.text = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_nested_wrap(subgroup~group_label, ncol = 6),
  theme(panel.spacing=unit(.075, "lines"))
) 

p5 + s
ggsave('plots/orc_spr_rt_subgroup3.png', width=6.5, height=2.25, dpi=600)

p6 + s
ggsave('plots/src_spr_rt_subgroup3.png', width=6.5, height=2.25, dpi=600)

# ---------------------------------------------------------------------------- #
# model: AJT (RP acceptance) vs. SPR (Log RT @ combined region) ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

temp <- ajt_ind2 %>%
  filter(subgroup %in% c('RP Rejector', 'Inconsistent'))

check <- temp %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

md <- spr_rt %>%
  filter(participant %in% temp$participant) %>%
  filter(region2 %in% c(1, 2, 3)) %>%
  group_by(study, group, participant, item, dependency, environment) %>%
  summarise(rt = sum(rt)) %>%
  ungroup() %>%
  select(study, group, participant, item, dependency, environment, rt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(rt = as.numeric(rt)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

check <- md %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, dependency, environment) %>%
  summarise(mean_rt = mean(rt, na.rm=T),
            lcl = my.lcl(rt),
            ucl = my.ucl(rt),
            sd = sd(rt, na.rm=T)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data = filter(plot, study == 'ORC'))
p2 <- ggplot(data = filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean_rt, group=dependency, col=dependency, linetype=dependency),
  geom_line(position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), linetype = 1, width = .5, position = position_dodge(width = .4)),
  geom_point(size = 1, shape = 21, fill = 'white', position = position_dodge(width = .4)),
  theme_classic(),
  scale_y_continuous(name="RT (ms)"),
  scale_x_discrete(name="Environment"),
  scale_colour_manual(name="Dependency", values=c('#648fff', '#ffb000'), labels=c('Gap', 'RP')),
  theme(
    text = element_text(size = 12),
    #legend.position = 'right',
    plot.margin = margin(2, 2, 2, 2),
    #legend.margin = margin(r = -.4, unit = 'cm')
  ),
  facet_manual(~group_label, design = 'ABC', axes = 'all', remove_labels = FALSE)
)

p1 + s

p2 + s

# table ...

table <- plot %>%
  select(study, Group = group, Environment = environment, 
         Dependency = dependency, Mean = mean_rt, LCL = lcl, UCL = ucl) %>%
  arrange(study, Group, Environment, Dependency) %>%
  mutate(Dependency = as.character(Dependency)) %>%
  mutate_if(is.character, str_replace_all, 
            c('gap' = 'Gap', 'pronoun' = 'RP'))

table %>% write_rds('tables/spr_rt_collapsed_nonacceptor.rds')

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model4 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + environment + group | item), 
               data = filter(md, study == 'ORC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model4)
toc()
beep()

# ... 170.41 sec elapsed

check_convergence(model4) # okay
check_singularity(model4) # warning (random effects table looks okay)
MuMIn::r.squaredGLMM(model4) # R2m = 0.0416322; R2c = 0.341798

# save ...

model4 %>% write_rds('models/orc_spr_logrt_subgroup_model4.rds')
model4 <- read_rds('models/orc_spr_logrt_subgroup_model4.rds')

# model summary ...

summary(model4)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)   
# (Intercept)                      6.930e+00  1.774e-02  1.399e+02 390.735  < 2e-16 ***
# dependency2                     -5.348e-02  1.103e-02  5.171e+01  -4.849 1.18e-05 ***
# environment2                    -1.253e-02  1.254e-02  6.652e+01  -1.000 0.320988    
# environment3                    -3.847e-03  1.290e-02  7.546e+01  -0.298 0.766354    
# group2                           4.942e-03  3.990e-02  1.382e+02   0.124 0.901609    
# group3                           1.376e-01  4.042e-02  1.358e+02   3.404 0.000872 ***
# dependency2:environment2        -8.631e-02  2.277e-02  3.562e+03  -3.790 0.000153 ***
# dependency2:environment3        -4.748e-02  2.273e-02  3.595e+03  -2.089 0.036745 *  
# dependency2:group2              -3.373e-02  2.343e-02  1.510e+02  -1.440 0.152043    
# dependency2:group3              -4.694e-02  2.394e-02  1.502e+02  -1.961 0.051706 .  
# environment2:group2              5.550e-02  2.698e-02  2.953e+02   2.057 0.040572 *  
# environment3:group2             -1.622e-03  2.789e-02  1.478e+02  -0.058 0.953713    
# environment2:group3             -1.803e-02  2.749e-02  2.922e+02  -0.656 0.512437    
# environment3:group3              2.564e-02  2.836e-02  1.448e+02   0.904 0.367310    
# dependency2:environment2:group2 -1.651e-01  5.185e-02  3.547e+03  -3.184 0.001465 ** 
# dependency2:environment3:group2  4.751e-02  5.158e-02  3.631e+03   0.921 0.357077    
# dependency2:environment2:group3 -1.193e-01  5.275e-02  3.627e+03  -2.261 0.023813 *  
# dependency2:environment3:group3  7.284e-02  5.270e-02  3.625e+03   1.382 0.167012  

# table ...

dt <- summary(model4)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

#dt %>% write_rds('tables/orc_spr_logrt_subgroup_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                              |6.93                     |0.02 |390.74 |< .001 |
# |Dependency[RP]                           |−0.05                    |0.01 |−4.85  |< .001 |
# |Environment[Long]                        |−0.01                    |0.01 |−1.00  |0.321  |
# |Environment[Island]                      |0.00                     |0.01 |−0.30  |0.766  |
# |Group[KLE]                               |0.00                     |0.04 |0.12   |0.902  |
# |Group[MLE]                               |0.14                     |0.04 |3.40   |< .001 |
# |Dependency[RP]:Environment[Long]         |−0.09                    |0.02 |−3.79  |< .001 |
# |Dependency[RP]:Environment[Island]       |−0.05                    |0.02 |−2.09  |0.037  |
# |Dependency[RP]:Group[KLE]                |−0.03                    |0.02 |−1.44  |0.152  |
# |Dependency[RP]:Group[MLE]                |−0.05                    |0.02 |−1.96  |0.052  |
# |Environment[Long]:Group[KLE]             |0.06                     |0.03 |2.06   |0.041  |
# |Environment[Island]:Group[KLE]           |0.00                     |0.03 |−0.06  |0.954  |
# |Environment[Long]:Group[MLE]             |−0.02                    |0.03 |−0.66  |0.512  |
# |Environment[Island]:Group[MLE]           |0.03                     |0.03 |0.90   |0.367  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |−0.17                    |0.05 |−3.18  |0.001  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.05                     |0.05 |0.92   |0.357  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |−0.12                    |0.05 |−2.26  |0.024  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.07                     |0.05 |1.38   |0.167  |

# joint tests ...

joint_tests(model4)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  23.509  <.0001 ***
# environment                    2 Inf   0.531  0.5881
# group                          2 Inf   6.219  0.0020 **
# dependency:environment         2 Inf   7.211  0.0007 **
# dependency:group               2 Inf   2.334  0.0969 .
# environment:group              4 Inf   2.861  0.0220 *
# dependency:environment:group   4 Inf   6.264  <.0001 ***

# post-hoc tests ...

pairwise <- model4 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group  estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS   -0.000207 0.0211 Inf  -0.010  1.0000
# pronoun - gap long        ENS    0.008261 0.0210 Inf   0.393  1.0000
# pronoun - gap island      ENS   -0.087809 0.0211 Inf  -4.171  0.0003 ***
# pronoun - gap short       KLE    0.005252 0.0321 Inf   0.164  1.0000
# pronoun - gap long        KLE   -0.151360 0.0320 Inf  -4.723  <.0001 ***
# pronoun - gap island      KLE   -0.034837 0.0322 Inf  -1.084  0.9461
# pronoun - gap short       MLE   -0.031674 0.0330 Inf  -0.961  0.9746
# pronoun - gap long        MLE   -0.142471 0.0329 Inf  -4.325  0.0001 ***
# pronoun - gap island      MLE   -0.046436 0.0330 Inf  -1.409  0.7872

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_subgroup_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:------|
# |ENS   |Short       |RP − Gap |0.00                     |0.02 |−0.01 |> .999 |
# |ENS   |Long        |RP − Gap |0.01                     |0.02 |0.39  |> .999 |
# |ENS   |Island      |RP − Gap |−0.09                    |0.02 |−4.17 |< .001 |
# |KLE   |Short       |RP − Gap |0.01                     |0.03 |0.16  |> .999 |
# |KLE   |Long        |RP − Gap |−0.15                    |0.03 |−4.72 |< .001 |
# |KLE   |Island      |RP − Gap |−0.03                    |0.03 |−1.08 |0.946  |
# |MLE   |Short       |RP − Gap |−0.03                    |0.03 |−0.96 |0.975  |
# |MLE   |Long        |RP − Gap |−0.14                    |0.03 |−4.32 |< .001 |
# |MLE   |Island      |RP − Gap |−0.05                    |0.03 |−1.41 |0.787  |

# plot ...

plot(emmeans(model4, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model8 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 + environment | item), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model8)
toc()
beep()

# ... 12.14 sec elapsed
# ... singular fit

# check performance ...

check_convergence(model8) # okay
check_singularity(model8) # warning (random effects table looks okay)
anova(model7, model8, refit = FALSE) # model 8 isn't better

# save ...

model8 %>% write_rds('models/src_spr_logrt_subgroup_model8.rds')
model8 <- read_rds('models/src_spr_logrt_subgroup_model8.rds')

# model summary ...

summary(model8)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                        7.19137    0.02828   93.77647 254.295  < 2e-16 ***
# dependency2                       -0.04658    0.01381  177.19214  -3.374 0.000909 ***
# environment2                       0.04041    0.01813   39.30124   2.230 0.031563 *  
# environment3                      -0.03883    0.01909   59.70730  -2.034 0.046393 *  
# group2                             0.17562    0.05485   90.68163   3.202 0.001885 ** 
# group3                             0.33191    0.07182   90.68583   4.621 1.26e-05 ***
# dependency2:environment2          -0.10347    0.03232 2486.90004  -3.202 0.001383 ** 
# dependency2:environment3          -0.11789    0.03225 2465.86966  -3.655 0.000262 ***
# dependency2:group2                -0.05487    0.02704  176.38639  -2.030 0.043885 *  
# dependency2:group3                -0.04565    0.03545  177.00702  -1.288 0.199505    
# environment2:group2                0.04636    0.03446   99.31061   1.345 0.181650    
# environment3:group2               -0.04516    0.03513  107.67276  -1.285 0.201412    
# environment2:group3                0.12437    0.04509   99.24429   2.758 0.006917 ** 
# environment3:group3                0.01199    0.04601  107.75612   0.261 0.794859    
# dependency2:environment2:group2   -0.18088    0.06331 2477.28238  -2.857 0.004311 ** 
# dependency2:environment3:group2    0.12497    0.06319 2461.33150   1.978 0.048063 *  
# dependency2:environment2:group3   -0.12987    0.08320 2493.72767  -1.561 0.118655    
# dependency2:environment3:group3    0.09349    0.08289 2466.18416   1.128 0.259467 

# table ...

dt <- summary(model8)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_subgroup_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                              |7.19                     |0.03 |254.30 |< .001 |
# |Dependency[RP]                           |−0.05                    |0.01 |−3.37  |< .001 |
# |Environment[Long]                        |0.04                     |0.02 |2.23   |0.032  |
# |Environment[Island]                      |−0.04                    |0.02 |−2.03  |0.046  |
# |Group[KLE]                               |0.18                     |0.05 |3.20   |0.002  |
# |Group[MLE]                               |0.33                     |0.07 |4.62   |< .001 |
# |Dependency[RP]:Environment[Long]         |−0.10                    |0.03 |−3.20  |0.001  |
# |Dependency[RP]:Environment[Island]       |−0.12                    |0.03 |−3.66  |< .001 |
# |Dependency[RP]:Group[KLE]                |−0.05                    |0.03 |−2.03  |0.044  |
# |Dependency[RP]:Group[MLE]                |−0.05                    |0.04 |−1.29  |0.200  |
# |Environment[Long]:Group[KLE]             |0.05                     |0.03 |1.35   |0.182  |
# |Environment[Island]:Group[KLE]           |−0.05                    |0.04 |−1.29  |0.201  |
# |Environment[Long]:Group[MLE]             |0.12                     |0.05 |2.76   |0.007  |
# |Environment[Island]:Group[MLE]           |0.01                     |0.05 |0.26   |0.795  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |−0.18                    |0.06 |−2.86  |0.004  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.12                     |0.06 |1.98   |0.048  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |−0.13                    |0.08 |−1.56  |0.119  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.09                     |0.08 |1.13   |0.259  |

# joint tests ...

joint_tests(model8)

# model term                   df1     df2 F.ratio p.value
# dependency                     1   90.68  11.346  0.0011 **
# environment                    2   38.02   8.296  0.0010 **
# group                          2   90.98  12.598  <.0001 ***
# dependency:environment         2 2381.41   7.911  0.0004 ***
# dependency:group               2   90.04   2.327  0.1034
# environment:group              4   90.31   3.128  0.0185 *
# dependency:environment:group   4 2373.18   6.362  <.0001 ***

# post-hoc tests ...

pairwise <- model8 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df t.ratio p.value
# pronoun - gap short       ENS     0.0505 0.0278 635   1.813  0.4794
# pronoun - gap long        ENS     0.0506 0.0279 636   1.816  0.4776
# pronoun - gap island      ENS    -0.1403 0.0279 636  -5.032  <.0001 ***
# pronoun - gap short       KLE     0.0142 0.0360 637   0.395  1.0000
# pronoun - gap long        KLE    -0.1665 0.0362 642  -4.606  <.0001 ***
# pronoun - gap island      KLE    -0.0515 0.0359 632  -1.434  0.7715
# pronoun - gap short       MLE     0.0169 0.0528 642   0.320  1.0000
# pronoun - gap long        MLE    -0.1128 0.0530 646  -2.128  0.2648
# pronoun - gap island      MLE    -0.0803 0.0527 635  -1.523  0.7072

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_subgroup_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$t$   |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:------|
# |ENS   |Short       |RP − Gap |0.05                     |0.03 |1.81  |0.479  |
# |ENS   |Long        |RP − Gap |0.05                     |0.03 |1.82  |0.478  |
# |ENS   |Island      |RP − Gap |−0.14                    |0.03 |−5.03 |< .001 |
# |KLE   |Short       |RP − Gap |0.01                     |0.04 |0.39  |> .999 |
# |KLE   |Long        |RP − Gap |−0.17                    |0.04 |−4.61 |< .001 |
# |KLE   |Island      |RP − Gap |−0.05                    |0.04 |−1.43 |0.771  |
# |MLE   |Short       |RP − Gap |0.02                     |0.05 |0.32  |> .999 |
# |MLE   |Long        |RP − Gap |−0.11                    |0.05 |−2.13 |0.265  |
# |MLE   |Island      |RP − Gap |−0.08                    |0.05 |−1.52 |0.707  |

# plot ...

plot(emmeans(model8, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... rrt src model (final)
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model7 <- lmer(rrt ~ dependency * environment * group + 
                 (1 + dependency | participant) + 
                 (1 + dependency | item), 
               data = filter(md, study == 'SRC'),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model7)
toc()
beep()

performance::check_convergence(model7) # TRUE
performance::check_singularity(model7) # FALSE

# 1.34 sec elapsed
# no warnings

# save ...

model7 %>% write_rds('models/src_spr_rrt_subgroup.rds')
model7 <- read_rds('models/src_spr_rrt_subgroup.rds')

# model summary ...

summary(model7)

# Parameter                       Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                       -56.78      14.99   64.03  -3.787 0.000339 ***
# dependency2                       -88.34      21.13   34.98  -4.181 0.000185 ***
# environment2                       68.56      24.17 2557.23   2.837 0.004595 ** 
# environment3                      119.08      24.16 2570.97   4.929 8.78e-07 ***
# group2                            107.68      27.69   91.49   3.889 0.000191 ***
# group3                            122.96      36.21   91.02   3.396 0.001015 ** 
# dependency2:environment2         -200.30      48.34 2556.31  -4.143 3.54e-05 ***
# dependency2:environment3         -161.89      48.35 2569.37  -3.349 0.000824 ***
# dependency2:group2                -81.95      40.63   91.63  -2.017 0.046596 *  
# dependency2:group3               -156.92      53.07   90.85  -2.957 0.003958 ** 
# environment2:group2                70.61      47.57 2578.40   1.484 0.137896    
# environment3:group2              -134.12      47.51 2578.49  -2.823 0.004794 ** 
# environment2:group3               165.23      62.06 2579.29   2.662 0.007809 ** 
# environment3:group3               -52.56      62.02 2564.47  -0.847 0.396798    
# dependency2:environment2:group2  -254.51      95.12 2579.15  -2.676 0.007502 ** 
# dependency2:environment3:group2    18.71      95.00 2579.34   0.197 0.843877    
# dependency2:environment2:group3  -248.51     124.24 2579.94  -2.000 0.045579 *  
# dependency2:environment3:group3   141.74     124.06 2564.15   1.143 0.253334 

# table ...

table <- summary(model7)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  rename('SE' = 'Std. Error', 't' = 't value', 'p' = 'Pr(>|t|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('dependency2' = 'Dependency[RP]',
              'environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'group2' = 'Group[KLE]', 
              'group3' = 'Group[MLE]'))

table %>% write_rds('tables/src_spr_rrt_subgroup_model.rds')

# joint tests ...

joint_tests(model7)

# model term                   df1     df2 F.ratio p.value
# dependency                     1   34.97  17.466  0.0002 ***
# environment                    2 2558.38  12.223  <.0001 ***
# group                          2   90.46  10.387  0.0001 ***
# dependency:environment         2 2557.46   9.658  0.0001 ***
# dependency:group               2   89.84   5.109  0.0079 **
# environment:group              4 2565.43   6.624  <.0001 ***
# dependency:environment:group   4 2565.11   4.037  0.0029 **

# post-hoc tests ...

pairwise <- model7 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate   SE  df t.ratio p.value
# pronoun - gap short       ENS       74.0 42.2 384   1.754  0.5263
# pronoun - gap long        ENS       41.3 42.2 384   0.980  0.9711
# pronoun - gap island      ENS     -141.4 42.0 379  -3.367  0.0075 **
# pronoun - gap short       KLE       70.6 54.2 519   1.303  0.8532
# pronoun - gap long        KLE     -216.5 54.2 520  -3.996  0.0007 ***
# pronoun - gap island      KLE     -126.1 54.1 517  -2.331  0.1669
# pronoun - gap short       MLE      -47.4 78.9 615  -0.600  0.9992
# pronoun - gap long        MLE     -328.5 78.9 615  -4.163  0.0003 ***
# pronoun - gap island      MLE     -121.0 78.9 615  -1.534  0.6993

# table ...

table <- pairwise %>%
  rename('Group' = 'group', 'Environment' = 'environment', 
         'Contrast' = 'contrast', 't' = 't.ratio', 'Estimate' = 'estimate',
         'p' = 'p.value') %>%
  select(Group, Environment, Contrast, everything()) %>%
  mutate_at(c('Environment', 'Contrast'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island',
              'gap' = 'Gap', 'pronoun' = 'RP')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Group, Environment)

table %>% write_rds('tables/src_spr_rrt_subgroup_posthoc.rds')

# plot ...

plot(emmeans(model7, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# model: AJT (RP acceptance) vs. SPR (Log RT @ single region) ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

temp <- ajt_ind2 %>%
  filter(pairing %in% c('EE', 'KE', 'ME')) %>%
  filter(subgroup %in% c('RP Rejector', 'Inconsistent'))

check <- temp %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

md <- spr_rt %>%
  filter(participant %in% temp$participant) %>%
  select(study, group, participant, item, dependency, environment, region2, rt) %>%
  mutate_at(c('study', 'group', 'participant', 'item', 'dependency', 'environment'), factor) %>%
  mutate(rt = as.numeric(rt)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

check <- md %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# plot ...

plot <- md %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')')) %>%
  group_by(study, group, group_label, dependency, environment, region2) %>%
  summarise(mean_rt = mean(log(rt), na.rm=T),
            lcl = my.lcl(log(rt)),
            ucl = my.ucl(log(rt)),
            sd = sd(rt, na.rm=T)) %>%
  ungroup() %>%
  mutate(environment = recode(environment, 'short' = 'Short', 'long' = 'Long', 'island' = 'Island'))

p1 <- ggplot(data = filter(plot, study == 'ORC', region2 == 1))
p2 <- ggplot(data = filter(plot, study == 'SRC', region2 == 2))

s <- list(
  aes(x=environment, y=mean_rt, group=dependency, col=dependency, 
      linetype=dependency),
  geom_line(position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                linetype = 1, width = .5, 
                position = position_dodge(width = .4)),
  geom_point(size = 1, shape = 21, fill = 'white', 
             position = position_dodge(width = .4),
             show.legend = FALSE),
  theme_classic(),
  scale_y_continuous(name="Log RT"),
  scale_x_discrete(name="Environment"),
  scale_colour_manual(name="Dependency", 
                      values=c('#648fff', '#ffb000'), 
                      labels=c('Gap', 'RP')),
  scale_linetype_manual(name="Dependency", 
                      values=c('42', 'solid'), 
                      labels=c('Gap', 'RP')),
  theme(
    text = element_text(size = 12),
    #legend.position = 'right',
    plot.margin = margin(2, 2, 2, 2),
    #legend.margin = margin(r = -.4, unit = 'cm')
  ),
  facet_manual(~group_label, design = 'ABC', axes = 'all', remove_labels = FALSE)
)

p1 + s + coord_cartesian(ylim = c(5.84, 6.25))

p2 + s + coord_cartesian(ylim = c(5.8, 6.52))

# table ...

table <- plot %>%
  select(study, Group = group, Environment = environment, 
         Dependency = dependency, Mean = mean_rt, LCL = lcl, UCL = ucl) %>%
  arrange(study, Group, Environment, Dependency) %>%
  mutate(Dependency = as.character(Dependency)) %>%
  mutate_if(is.character, str_replace_all, 
            c('gap' = 'Gap', 'pronoun' = 'RP'))

table %>% write_rds('tables/spr_rt_r2_nonacceptor.rds')

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model6 <- lmer(log(rt) ~ dependency * environment * group + 
                 (1 + dependency + environment | participant) + 
                 (1 + dependency + group | item), 
               data = filter(md, study == 'ORC', region2 == 1),
               control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model6)
toc()
beep()

# ... 30.27 sec elapsed
# ... singular fit

check_convergence(model6) # okay
check_singularity(model6) # warning (random effects table looks okay)
MuMIn::r.squaredGLMM(model6) # R2m = 0.04068912; R2c = 0.3290754

# save ...

model6 %>% write_rds('models/orc_spr_logrt_subgroup_model.rds')
model6 <- read_rds('models/orc_spr_logrt_subgroup_model.rds')

# model summary ...

summary(model6)

# Parameter                       Estimate Std. Error       df t value Pr(>|t|)   
# (Intercept)                      6.000e+00  2.102e-02  1.397e+02 285.496  < 2e-16 ***
# dependency2                     -7.106e-02  1.479e-02  5.042e+01  -4.804 1.44e-05 ***
# environment2                    -6.421e-03  1.563e-02  1.321e+02  -0.411 0.681946    
# environment3                    -1.642e-02  1.471e-02  1.303e+02  -1.116 0.266433    
# group2                           6.156e-02  4.609e-02  1.355e+02   1.336 0.183901    
# group3                           1.719e-01  4.798e-02  1.369e+02   3.582 0.000472 ***
# dependency2:environment2        -9.115e-02  2.855e-02  3.374e+03  -3.193 0.001420 ** 
# dependency2:environment3        -1.693e-02  2.834e-02  3.353e+03  -0.597 0.550274    
# dependency2:group2              -3.055e-02  3.053e-02  1.312e+02  -1.001 0.318826    
# dependency2:group3              -4.589e-02  3.102e-02  1.287e+02  -1.479 0.141519    
# environment2:group2              2.616e-02  3.555e-02  1.322e+02   0.736 0.463195    
# environment3:group2             -1.326e-02  3.354e-02  1.315e+02  -0.395 0.693230    
# environment2:group3             -3.674e-02  3.604e-02  1.294e+02  -1.020 0.309830    
# environment3:group3              2.772e-04  3.398e-02  1.286e+02   0.008 0.993504    
# dependency2:environment2:group2 -9.628e-02  6.475e-02  3.311e+03  -1.487 0.137111    
# dependency2:environment3:group2  6.937e-02  6.452e-02  3.376e+03   1.075 0.282397    
# dependency2:environment2:group3 -1.688e-01  6.573e-02  3.383e+03  -2.568 0.010272 *  
# dependency2:environment3:group3  1.226e-01  6.554e-02  3.354e+03   1.871 0.061394 . 

# table ...

dt <- summary(model6)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_subgroup_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_subgroup_model.csv')
dt <- read_rds('tables/orc_spr_logrt_subgroup_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                              |6.00                     |0.02 |285.50 |< .001 |
# |Dependency[RP]                           |−0.07                    |0.01 |−4.80  |< .001 |
# |Environment[Long]                        |−0.01                    |0.02 |−0.41  |0.682  |
# |Environment[Island]                      |−0.02                    |0.01 |−1.12  |0.266  |
# |Group[KLE]                               |0.06                     |0.05 |1.34   |0.184  |
# |Group[MLE]                               |0.17                     |0.05 |3.58   |< .001 |
# |Dependency[RP]:Environment[Long]         |−0.09                    |0.03 |−3.19  |0.001  |
# |Dependency[RP]:Environment[Island]       |−0.02                    |0.03 |−0.60  |0.550  |
# |Dependency[RP]:Group[KLE]                |−0.03                    |0.03 |−1.00  |0.319  |
# |Dependency[RP]:Group[MLE]                |−0.05                    |0.03 |−1.48  |0.142  |
# |Environment[Long]:Group[KLE]             |0.03                     |0.04 |0.74   |0.463  |
# |Environment[Island]:Group[KLE]           |−0.01                    |0.03 |−0.40  |0.693  |
# |Environment[Long]:Group[MLE]             |−0.04                    |0.04 |−1.02  |0.310  |
# |Environment[Island]:Group[MLE]           |0.00                     |0.03 |0.01   |0.994  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |−0.10                    |0.06 |−1.49  |0.137  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.07                     |0.06 |1.08   |0.282  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |−0.17                    |0.07 |−2.57  |0.010  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.12                     |0.07 |1.87   |0.061  |

# joint tests ...

joint_tests(model6)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  23.080  <.0001 ***
# environment                    2 Inf   0.638  0.5284
# group                          2 Inf   6.468  0.0016 **
# dependency:environment         2 Inf   5.755  0.0032 **
# dependency:group               2 Inf   1.274  0.2797
# environment:group              4 Inf   0.941  0.4386
# dependency:environment:group   4 Inf   5.447  0.0002 ***

# post-hoc tests ...

pairwise <- model6 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE  df z.ratio p.value
# pronoun - gap short       ENS   -0.01767 0.0269 Inf  -0.657  0.9983
# pronoun - gap long        ENS   -0.02045 0.0267 Inf  -0.765  0.9947
# pronoun - gap island      ENS   -0.09860 0.0268 Inf  -3.684  0.0021 **
# pronoun - gap short       KLE   -0.03925 0.0410 Inf  -0.957  0.9747
# pronoun - gap long        KLE   -0.13832 0.0412 Inf  -3.354  0.0071 **
# pronoun - gap island      KLE   -0.05082 0.0411 Inf  -1.236  0.8854
# pronoun - gap short       MLE   -0.04817 0.0420 Inf  -1.146  0.9242
# pronoun - gap long        MLE   -0.21975 0.0419 Inf  -5.248  <.0001 ***
# pronoun - gap island      MLE   -0.00647 0.0417 Inf  -0.155  1.0000

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_logrt_subgroup_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_logrt_subgroup_posthoc.csv')
dt <- read_rds('tables/orc_spr_logrt_subgroup_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:------|
# |ENS   |Short       |RP − Gap |−0.02                    |0.03 |−0.66 |0.998  |
# |ENS   |Long        |RP − Gap |−0.02                    |0.03 |−0.76 |0.995  |
# |ENS   |Island      |RP − Gap |−0.10                    |0.03 |−3.68 |0.002  |
# |KLE   |Short       |RP − Gap |−0.04                    |0.04 |−0.96 |0.975  |
# |KLE   |Long        |RP − Gap |−0.14                    |0.04 |−3.35 |0.007  |
# |KLE   |Island      |RP − Gap |−0.05                    |0.04 |−1.24 |0.885  |
# |MLE   |Short       |RP − Gap |−0.05                    |0.04 |−1.15 |0.924  |
# |MLE   |Long        |RP − Gap |−0.22                    |0.04 |−5.25 |< .001 |
# |MLE   |Island      |RP − Gap |−0.01                    |0.04 |−0.16 |> .999 |

# plot ...

plot(emmeans(model6, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model11 <- lmer(log(rt) ~ dependency * environment * group + 
                  (1 + dependency * environment | participant) + 
                  (1 + dependency + environment | item), 
                data = filter(md, study == 'SRC', region2 == 2),
                control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e6)))
summary(model11)
toc()
beep()

# ... 664.93 sec elapsed
# ... singular fit

check_convergence(model11) # warning (lme4 says it's okay)
check_singularity(model11) # warning (random effects table looks okay)

# save ...

model11 %>% write_rds('models/src_spr_logrt_subgroup_model.rds')
model11 <- read_rds('models/src_spr_logrt_subgroup_model.rds')

# model summary ...

summary(model11)

# Parameter                         Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                       6.07246    0.03219  88.33274 188.631  < 2e-16 ***
# dependency2                      -0.05984    0.01561 119.87824  -3.832 0.000204 ***
# environment2                      0.04478    0.02030  93.06146   2.206 0.029880 *  
# environment3                     -0.03910    0.02062  55.14831  -1.897 0.063106 .  
# group2                            0.16164    0.06230  84.89013   2.595 0.011150 *  
# group3                            0.25568    0.08150  84.66997   3.137 0.002348 ** 
# dependency2:environment2         -0.18861    0.04145 104.62008  -4.551 1.45e-05 ***
# dependency2:environment3         -0.10021    0.03801 197.60509  -2.637 0.009040 ** 
# dependency2:group2               -0.02716    0.02977 765.06884  -0.912 0.361883    
# dependency2:group3               -0.01156    0.03875 754.48266  -0.298 0.765584    
# environment2:group2               0.06645    0.03987 111.51516   1.667 0.098352 .  
# environment3:group2              -0.07194    0.03862 153.60161  -1.863 0.064392 .  
# environment2:group3               0.17462    0.05161 107.79388   3.383 0.000999 ***
# environment3:group3              -0.01961    0.05041 152.73758  -0.389 0.697882    
# dependency2:environment2:group2  -0.15428    0.08172 105.84786  -1.888 0.061763 .  
# dependency2:environment3:group2   0.22543    0.07454 196.31454   3.024 0.002827 ** 
# dependency2:environment2:group3  -0.13943    0.10617 103.35564  -1.313 0.191988    
# dependency2:environment3:group3   0.25852    0.09742 196.32586   2.654 0.008611 ** 

# table ...

dt <- summary(model11)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_subgroup_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_subgroup_model.csv')
dt <- read_rds('tables/src_spr_logrt_subgroup_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                |Estimate~($\hat{\beta}$) |$SE$ |$t$    |$p$    |
# |:----------------------------------------|:------------------------|:----|:------|:------|
# |(Intercept)                              |6.07                     |0.03 |188.63 |< .001 |
# |Dependency[RP]                           |−0.06                    |0.02 |−3.83  |< .001 |
# |Environment[Long]                        |0.04                     |0.02 |2.21   |0.030  |
# |Environment[Island]                      |−0.04                    |0.02 |−1.90  |0.063  |
# |Group[KLE]                               |0.16                     |0.06 |2.59   |0.011  |
# |Group[MLE]                               |0.26                     |0.08 |3.14   |0.002  |
# |Dependency[RP]:Environment[Long]         |−0.19                    |0.04 |−4.55  |< .001 |
# |Dependency[RP]:Environment[Island]       |−0.10                    |0.04 |−2.64  |0.009  |
# |Dependency[RP]:Group[KLE]                |−0.03                    |0.03 |−0.91  |0.362  |
# |Dependency[RP]:Group[MLE]                |−0.01                    |0.04 |−0.30  |0.766  |
# |Environment[Long]:Group[KLE]             |0.07                     |0.04 |1.67   |0.098  |
# |Environment[Island]:Group[KLE]           |−0.07                    |0.04 |−1.86  |0.064  |
# |Environment[Long]:Group[MLE]             |0.17                     |0.05 |3.38   |< .001 |
# |Environment[Island]:Group[MLE]           |−0.02                    |0.05 |−0.39  |0.698  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |−0.15                    |0.08 |−1.89  |0.062  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |0.23                     |0.07 |3.02   |0.003  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |−0.14                    |0.11 |−1.31  |0.192  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.26                     |0.10 |2.65   |0.009  |

# joint tests ...

joint_tests(model11)

# model term                   df1   df2 F.ratio p.value
# dependency                     1 31.77  14.610  0.0006 ***
# environment                    2 37.63   6.494  0.0038 **
# group                          2 90.71   6.483  0.0023 **
# dependency:environment         2 90.23  10.406  0.0001 ***
# dependency:group               2 89.58   0.415  0.6613
# environment:group              4 88.75   5.253  0.0008 ***
# dependency:environment:group   4 89.15   7.853  <.0001 ***

# post-hoc tests ...

pairwise <- model11 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate     SE   df t.ratio p.value
# pronoun - gap short       ENS     0.0705 0.0328 80.7   2.150  0.2655
# pronoun - gap long        ENS    -0.0202 0.0343 82.3  -0.589  0.9992
# pronoun - gap island      ENS    -0.1910 0.0324 81.7  -5.896  <.0001 ***
# pronoun - gap short       KLE     0.0196 0.0433 92.7   0.453  0.9999
# pronoun - gap long        KLE    -0.2254 0.0452 92.5  -4.983  <.0001 ***
# pronoun - gap island      KLE    -0.0165 0.0418 88.6  -0.395  1.0000
# pronoun - gap short       MLE     0.0192 0.0623 91.1   0.309  1.0000
# pronoun - gap long        MLE    -0.2109 0.0648 88.7  -3.253  0.0144 *
# pronoun - gap island      MLE     0.0162 0.0613 92.4   0.265  1.0000

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_logrt_subgroup_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_logrt_subgroup_posthoc.csv')
dt <- read_rds('tables/src_spr_logrt_subgroup_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$t$   |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:-----|:------|
# |ENS   |Short       |RP − Gap |0.07                     |0.03 |2.15  |0.265  |
# |ENS   |Long        |RP − Gap |−0.02                    |0.03 |−0.59 |> .999 |
# |ENS   |Island      |RP − Gap |−0.19                    |0.03 |−5.90 |< .001 |
# |KLE   |Short       |RP − Gap |0.02                     |0.04 |0.45  |> .999 |
# |KLE   |Long        |RP − Gap |−0.23                    |0.05 |−4.98 |< .001 |
# |KLE   |Island      |RP − Gap |−0.02                    |0.04 |−0.40 |> .999 |
# |MLE   |Short       |RP − Gap |0.02                     |0.06 |0.31  |> .999 |
# |MLE   |Long        |RP − Gap |−0.21                    |0.06 |−3.25 |0.014  |
# |MLE   |Island      |RP − Gap |0.02                     |0.06 |0.26  |> .999 |

# plot ...

plot(emmeans(model11, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: AJT (RP acceptance) vs. SPR (accuracy) ####
# ---------------------------------------------------------------------------- #

ajt_ind2 <- ajt_ind %>%
  group_by(study, group, language, pairing, participant, proficiency) %>%
  summarise(max = max(acceptance)) %>%
  ungroup() %>%
  mutate(subgroup = case_when(max >= .8 ~ "RP Acceptor",
                              max <= .2 ~ "RP Rejector",
                              TRUE ~ "Inconsistent")) %>%
  filter(language == 'English')

check <- ajt_ind2 %>%
  group_by(study, pairing) %>%
  summarise(n = n()) %>%
  ungroup()

sp2 <- spr_clean %>%
  filter(dependency %in% c('gap', 'pronoun')) %>%
  mutate(accuracy = as.logical(accuracy)) %>%
  group_by(study, group, participant, item, dependency, environment, accuracy) %>%
  summarise() %>%
  ungroup() %>%
  left_join(ajt_ind2, by = c('study', 'group', 'participant')) %>%
  filter(subgroup %in% c('RP Acceptor', 'RP Rejector', 'Inconsistent')) %>%
  mutate(subgroup = fct_relevel(subgroup, 'RP Acceptor', 'RP Rejector', 'Inconsistent')) %>%
  mutate(subgroup2 = case_when(subgroup %in% c('RP Rejector', 'Inconsistent') ~ 'RP Non-Acceptor',
                               subgroup == 'RP Acceptor' ~ 'RP Acceptor')) %>%
  mutate(subgroup2 = fct_relevel(subgroup2, 'RP Acceptor', 'RP Non-Acceptor'))

plot <- sp2 %>%
  group_by(study, group, subgroup) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, ' (*n* = ', group_size, ')')) %>%
  group_by(study, group, group_label, subgroup, dependency, environment) %>%
  summarise(mean = mean(accuracy, na.rm=T) * 100,
            lcl = my.lcl(accuracy) * 100,
            ucl = my.ucl(accuracy) * 100) %>%
  ungroup() %>%
  mutate(lcl = case_when(is.na(lcl) ~ mean, TRUE ~ lcl),
         ucl = case_when(is.na(ucl) ~ mean, TRUE ~ ucl))

plot2 <- sp2 %>%
  mutate(subgroup = subgroup2) %>%
  group_by(study, group, subgroup) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, ' (*n* = ', group_size, ')')) %>%
  group_by(study, group, group_label, subgroup, dependency, environment) %>%
  summarise(mean = mean(accuracy, na.rm=T) * 100,
            lcl = my.lcl(accuracy) * 100,
            ucl = my.ucl(accuracy) * 100) %>%
  ungroup() %>%
  mutate(lcl = case_when(is.na(lcl) ~ mean, TRUE ~ lcl),
         ucl = case_when(is.na(ucl) ~ mean, TRUE ~ ucl)) %>%
  mutate(subgroup2 = case_when(subgroup == 'RP Acceptor' ~ '+ RP Acceptor',
                               TRUE ~ '\u2212 RP Acceptor'))

# define data for plots ...

p3a <- ggplot(data=filter(plot2, study == 'ORC', subgroup == 'RP Acceptor'))
p3b <- ggplot(data=filter(plot2, study == 'ORC', subgroup == 'RP Non-Acceptor'))
p4a <- ggplot(data=filter(plot2, study == 'SRC', subgroup == 'RP Acceptor'))
p4b <- ggplot(data=filter(plot2, study == 'SRC', subgroup == 'RP Non-Acceptor'))

# styling for just non-acceptors ...

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency),
            position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1,
             shape = 21,
             fill = 'white',
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(50, 100)),
  scale_y_continuous(name='Tokens',
                     breaks = c(50, 75, 100)),
  scale_x_discrete(name='Environment',
                   limits = c('short', 'long', 'island'),
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency',
                      values=c('#648fff', '#ffb000'),
                      breaks = c('gap', 'pronoun'),
                      labels = c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency',
                        values=c('41', 'solid'),
                        breaks = c('gap', 'pronoun'),
                        labels = c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'right',
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

p3b + s
ggsave('plots/orc_spr_accuracy_nonacceptor.png', width=6.5, height=1.75, dpi=600)

p4b + s
ggsave('plots/src_spr_accuracy_nonacceptor.png', width=6.5, height=1.75, dpi=600)

# styling ...

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency),
            position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1,
             shape = 21,
             fill = 'white',
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(50, 100)),
  scale_y_continuous(name='Tokens',
                     breaks = c(50, 75, 100)),
  scale_x_discrete(name='Environment',
                   limits = c('short', 'long', 'island'),
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency',
                      values=c('#648fff', '#ffb000'),
                      breaks = c('gap', 'pronoun'),
                      labels = c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency',
                     values=c('41', 'solid'),
                     breaks = c('gap', 'pronoun'),
                     labels = c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        strip.text = element_markdown(),
        legend.position = 'right',
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_grid2(subgroup2~group_label,
              axes = 'all',
              remove_labels = FALSE)
)

plot3a = p3a + s + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
plot3b = p3b + s + theme(axis.title.y = element_blank())
grid.arrange(patchworkGrob(plot3a / plot3b + 
                             plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                             plot_layout(guides = 'collect')), left = 'Accuracy (%)')
g <- arrangeGrob(patchworkGrob(plot3a / plot3b + 
                                 plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                                 plot_layout(guides = 'collect')), left = 'Accuracy (%)')
ggsave('plots/orc_spr_accuracy_subgroup2.png', g, width=6.5, height=3.25, dpi=600)

plot4a = p4a + s + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
plot4b = p4b + s + theme(axis.title.y = element_blank())
grid.arrange(patchworkGrob(plot4a / plot4b + 
                             plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                             plot_layout(guides = 'collect')), left = 'Accuracy (%)')
g <- arrangeGrob(patchworkGrob(plot4a / plot4b + 
                                 plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))) +
                                 plot_layout(guides = 'collect')), left = 'Accuracy (%)')
ggsave('plots/src_spr_accuracy_subgroup2.png', g, width=6.5, height=3.25, dpi=600)

# nested plots ...

p5 <- ggplot(data=filter(plot2, study == 'ORC'))
p6 <- ggplot(data=filter(plot2, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=dependency, col=dependency, shape=dependency),
  geom_line(aes(linetype = dependency), position = position_dodge(width = .6)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
                linetype = 1,
                width = .5,
                position = position_dodge(width = .6),
                show.legend = FALSE),
  geom_point(size = 1,
             shape = 21,
             fill = 'white',
             position = position_dodge(width = .6),
             show.legend = FALSE),
  theme_classic(),
  coord_cartesian(ylim = c(49, 101)),
  scale_y_continuous(name='Accuracy (%)', breaks = c(50, 75, 100)),
  scale_x_discrete(name='Environment',
                   labels = c('Short', 'Long', 'Island')),
  scale_colour_manual(name='Dependency:',
                      values=c('#648fff', '#ffb000'),
                      breaks = c('gap', 'pronoun'),
                      labels = c('Gap', 'RP')),
  scale_linetype_manual(name='Dependency:',
                        values=c('41', 'solid'),
                        breaks = c('gap', 'pronoun'),
                        labels = c('Gap', 'RP')),
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom',
        strip.text = element_markdown(),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(t = -.4, unit = 'cm')),
  facet_nested_wrap(subgroup~group_label, ncol = 6),
  theme(panel.spacing=unit(.075, "lines"))
) 

p5 + s
ggsave('plots/orc_spr_accuracy_subgroup3.png', width=6.5, height=2.25, dpi=600)

p6 + s
ggsave('plots/src_spr_accuracy_subgroup3.png', width=6.5, height=2.25, dpi=600)

# ---------------------------------------------------------------------------- #
# model: AJT (RP acceptance) vs. SPR (accuracy) ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

temp <- ajt_ind2 %>%
  filter(subgroup %in% c('RP Rejector', 'Inconsistent'))

check <- temp %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant))) %>%
  ungroup()

md <- spr_clean %>%
  filter(participant %in% temp$participant) %>%
  filter(class == 'critical') %>%
  group_by(study, group, participant) %>%
  mutate(acc_rate = mean(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(study, group, participant, item, dependency, environment, acc_rate, accuracy) %>%
  summarise() %>%
  ungroup() %>%
  mutate(dependency = fct_drop(dependency),
         environment = fct_drop(environment),
         accuracy = as.logical(accuracy)) %>%
  mutate_at(c('study', 'group', 'participant', 'item'), factor) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

check <- md %>%
  group_by(study, group) %>%
  summarise(group_size = length(unique(participant)))

# apply deviation coding ...

contrasts(md$dependency) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

# view contrasts ...

contrasts(md$dependency)
contrasts(md$environment)
contrasts(md$group)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model15b <- glmer(accuracy ~ dependency * environment * group + 
                    (1 + dependency * environment | participant) + 
                    (1 + dependency + environment + group | item), 
                  data = filter(md, study == 'ORC'), family = binomial, 
                  control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model15b)
toc()
beep()

# ... 4885.7 sec elapsed

# check performance ...

check_convergence(model15b) # okay
check_singularity(model15b) # warning (random effects table looks okay)

# save ...

model15b %>% write_rds('models/orc_spr_accuracy_subgroup_model15b.rds')
model15b <- read_rds('models/orc_spr_accuracy_subgroup_model15b.rds')

# table ...

dt <- summary(model15b)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/orc_spr_accuracy_subgroup_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_accuracy_subgroup_model.csv')
dt <- read_rds('tables/orc_spr_accuracy_subgroup_model.rds')

dt %>% kbl(format = 'markdown')

# |Parameter                                     | Estimate|   SE|     z|     p|
# |:---------------------------------------------|--------:|----:|-----:|-----:|
# |(Intercept)                                   |     3.41| 0.21| 16.06| 0.000| ***
# |Dependency[RP]                                |     1.01| 0.32|  3.20| 0.001| **
# |Environment[Long]                             |     0.18| 0.39|  0.47| 0.640|
# |Environment[Island]                           |     0.21| 0.37|  0.57| 0.569|
# |Group[KLE]                                    |    -1.27| 0.34| -3.67| 0.000| ***
# |Group[MLE]                                    |    -0.97| 0.37| -2.64| 0.008| **
# |Dependency[RP]:Environment[Long]              |     0.70| 0.71|  0.99| 0.321|
# |Dependency[RP]:Environment[Island]            |     0.72| 0.73|  0.99| 0.324|
# |Dependency[RP]:Group[KLE]                     |     0.98| 0.40|  2.45| 0.014| *
# |Dependency[RP]:Group[MLE]                     |    -0.07| 0.40| -0.18| 0.855|
# |Environment[Long]:Group[KLE]                  |     0.18| 0.50|  0.37| 0.713|
# |Environment[Island]:Group[KLE]                |     0.67| 0.51|  1.30| 0.193|
# |Environment[Long]:Group[MLE]                  |     0.24| 0.51|  0.46| 0.643|
# |Environment[Island]:Group[MLE]                |     0.65| 0.52|  1.25| 0.211|
# |Dependency[RP]:Environment[Long]:Group[KLE]   |    -0.18| 0.97| -0.18| 0.855|
# |Dependency[RP]:Environment[Island]:Group[KLE] |     0.04| 0.94|  0.04| 0.968|
# |Dependency[RP]:Environment[Long]:Group[MLE]   |    -0.54| 0.98| -0.55| 0.581|
# |Dependency[RP]:Environment[Island]:Group[MLE] |    -0.40| 0.94| -0.42| 0.672|

# joint tests ...

joint_tests(model15b)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  10.262  0.0014 **
# environment                    2 Inf   0.198  0.8202
# group                          2 Inf   7.767  0.0004 ***
# dependency:environment         2 Inf   0.718  0.4879
# dependency:group               2 Inf   3.940  0.0194 *
# environment:group              4 Inf   0.608  0.6565
# dependency:environment:group   4 Inf   0.102  0.9820

# post-hoc tests ...

pairwise <- model15b %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/orc_spr_accuracy_subgroup_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/orc_spr_accuracy_subgroup_posthoc.csv')
dt <- read_rds('tables/orc_spr_accuracy_subgroup_posthoc.rds')

dt %>% kbl(format = 'markdown')

# |Group |Environment |Contrast |Estimate~($\hat{\beta}$) |$SE$ |$z$  |$p$    |
# |:-----|:-----------|:--------|:------------------------|:----|:----|:------|
# |ENS   |Short       |RP − Gap |0.12                     |0.62 |0.19 |> .999 |
# |ENS   |Long        |RP − Gap |1.06                     |0.70 |1.51 |0.665  |
# |ENS   |Island      |RP − Gap |0.96                     |0.68 |1.41 |0.740  |
# |KLE   |Short       |RP − Gap |1.14                     |0.54 |2.13 |0.241  |
# |KLE   |Long        |RP − Gap |1.90                     |0.65 |2.94 |0.028  | *
# |KLE   |Island      |RP − Gap |2.02                     |0.75 |2.69 |0.059  | .
# |MLE   |Short       |RP − Gap |0.36                     |0.58 |0.62 |0.998  |
# |MLE   |Long        |RP − Gap |0.76                     |0.75 |1.01 |0.949  |
# |MLE   |Island      |RP − Gap |0.80                     |0.74 |1.08 |0.926  |

# plot ...

plot(emmeans(model15b, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

# fit model ...

tic()
model9 <- glmer(accuracy ~ dependency * environment * group + 
                  (1 + dependency * environment | participant) + 
                  (1 + dependency * environment + group | item), 
                data = filter(md, study == 'SRC'), family = binomial, 
                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6)))
summary(model9)
toc()
beep()

# ... 6394.46 sec elapsed (107 min)

check_convergence(model9) # okay
check_singularity(model9) # warning (but random effects table looks okay)

# save ...

model9 %>% write_rds('models/src_spr_accuracy_subgroup_model9.rds')
model9 <- read_rds('models/src_spr_accuracy_subgroup_model9.rds')

# model summary ...

summary(model9)

# Parameter                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                       2.9937     0.2382  12.567  < 2e-16 ***
# dependency2                       1.6562     0.3380   4.900 9.58e-07 ***
# environment2                     -0.6191     0.4796  -1.291   0.1967    
# environment3                     -0.9878     0.4394  -2.248   0.0246 *  
# group2                           -0.2993     0.3516  -0.851   0.3946    
# group3                           -0.6623     0.4620  -1.434   0.1517    
# dependency2:environment2          0.2656     0.9152   0.290   0.7717    
# dependency2:environment3          0.3567     0.8655   0.412   0.6803    
# dependency2:group2                0.3406     0.4005   0.850   0.3952    
# dependency2:group3               -0.2084     0.4783  -0.436   0.6630    
# environment2:group2              -0.5998     0.5305  -1.131   0.2582    
# environment3:group2              -0.7742     0.4898  -1.581   0.1140    
# environment2:group3               0.2740     0.5985   0.458   0.6470    
# environment3:group3               0.8746     0.5794   1.509   0.1312    
# dependency2:environment2:group2  -0.0336     0.9879  -0.034   0.9729    
# dependency2:environment3:group2  -1.2244     0.9188  -1.333   0.1827    
# dependency2:environment2:group3   0.7232     1.1206   0.645   0.5187    
# dependency2:environment3:group3   0.2798     1.1202   0.250   0.8027 

# table ...

dt <- summary(model9)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  my_table()

dt %>% write_rds('tables/src_spr_accuracy_subgroup_model.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_accuracy_subgroup_model.csv')
dt <- read_rds('tables/src_spr_accuracy_subgroup_model.rds')

table %>% kbl(format = 'markdown')

# |Parameter                                |Estimate~($\hat{\beta}$) |$SE$ |$z$   |$p$    |
# |:----------------------------------------|:------------------------|:----|:-----|:------|
# |(Intercept)                              |2.99                     |0.24 |12.57 |< .001 |
# |Dependency[RP]                           |1.66                     |0.34 |4.90  |< .001 |
# |Environment[Long]                        |−0.62                    |0.48 |−1.29 |0.197  |
# |Environment[Island]                      |−0.99                    |0.44 |−2.25 |0.025  |
# |Group[KLE]                               |−0.30                    |0.35 |−0.85 |0.395  |
# |Group[MLE]                               |−0.66                    |0.46 |−1.43 |0.152  |
# |Dependency[RP]:Environment[Long]         |0.27                     |0.92 |0.29  |0.772  |
# |Dependency[RP]:Environment[Island]       |0.36                     |0.87 |0.41  |0.680  |
# |Dependency[RP]:Group[KLE]                |0.34                     |0.40 |0.85  |0.395  |
# |Dependency[RP]:Group[MLE]                |−0.21                    |0.48 |−0.44 |0.663  |
# |Environment[Long]:Group[KLE]             |−0.60                    |0.53 |−1.13 |0.258  |
# |Environment[Island]:Group[KLE]           |−0.77                    |0.49 |−1.58 |0.114  |
# |Environment[Long]:Group[MLE]             |0.27                     |0.60 |0.46  |0.647  |
# |Environment[Island]:Group[MLE]           |0.87                     |0.58 |1.51  |0.131  |
# |Dependency[RP]:Envir.[Long]:Group[KLE]   |−0.03                    |0.99 |−0.03 |0.973  |
# |Dependency[RP]:Envir.[Island]:Group[KLE] |−1.22                    |0.92 |−1.33 |0.183  |
# |Dependency[RP]:Envir.[Long]:Group[MLE]   |0.72                     |1.12 |0.65  |0.519  |
# |Dependency[RP]:Envir.[Island]:Group[MLE] |0.28                     |1.12 |0.25  |0.803  |

# joint tests ...

joint_tests(model9)

# model term                   df1 df2 F.ratio p.value
# dependency                     1 Inf  24.012  <.0001
# environment                    2 Inf   2.566  0.0768
# group                          2 Inf   1.110  0.3295
# dependency:environment         2 Inf   0.085  0.9182
# dependency:group               2 Inf   0.633  0.5312
# environment:group              4 Inf   1.775  0.1306
# dependency:environment:group   4 Inf   0.778  0.5393

# post-hoc tests ...

pairwise <- model9 %>%
  emmeans(~ dependency * environment * group) %>%
  contrast('revpairwise', simple = 'dependency', combine = TRUE) %>%
  summary(by = NULL, adjust = 'mvt')
pairwise

# contrast      environment group estimate    SE  df z.ratio p.value
# pronoun - gap short       ENS      1.376 0.792 Inf   1.738  0.4949
# pronoun - gap long        ENS      1.412 0.650 Inf   2.173  0.2187
# pronoun - gap island      ENS      2.048 0.550 Inf   3.725  0.0017 **
# pronoun - gap short       KLE      2.136 0.983 Inf   2.172  0.2191
# pronoun - gap long        KLE      2.138 0.700 Inf   3.053  0.0193 *
# pronoun - gap island      KLE      1.583 0.556 Inf   2.849  0.0367 *
# pronoun - gap short       MLE      0.834 0.812 Inf   1.027  0.9435
# pronoun - gap long        MLE      1.593 0.785 Inf   2.030  0.2955
# pronoun - gap island      MLE      1.785 0.821 Inf   2.175  0.2178

# table ...

dt <- pairwise %>%
  select(group, environment, contrast, everything()) %>%
  select(-df) %>%
  my_table()

dt %>% write_rds('tables/src_spr_accuracy_subgroup_posthoc.rds')
csv <- dt %>% csv_table() %>% write_csv('tables/src_spr_accuracy_subgroup_posthoc.csv')
dt <- read_rds('tables/src_spr_accuracy_subgroup_posthoc.rds')

# plot ...

plot(emmeans(model9, ~ dependency * environment | group), 
     combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

# ---------------------------------------------------------------------------- #
# plot: AJT (RP rating) vs. EPT (token) ####
# ---------------------------------------------------------------------------- #

# prep dataframe for plotting ...

aj2 <- ajt_clean %>%
  filter(language == 'English', 
         dependency == 'pronoun',
         environment %in% c('short', 'long', 'island')) %>%
  group_by(study, group, participant, environment) %>%
  summarise(rating = mean(response, na.rm = TRUE),
            zscore = mean(zscore, na.rm = TRUE)) %>%
  ungroup()

ep2 <- ep %>%
  filter(is.na(condition) == FALSE) %>%
  group_by(study, group, participant, environment) %>%
  summarise(token = length(which(type == 'resumption'))) %>%
  ungroup()

plot <- ep2 %>%
  left_join(aj2, by = c('study', 'group', 'participant', 'environment'))

plot <- plot %>%
  filter(is.na(rating) == FALSE) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island')) %>%
  mutate(environment = fct_recode(environment, 'Short' = 'short', 'Long' = 'long', 'Island' = 'island')) %>%
  group_by(study, group) %>%
  mutate(group_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(group_label = paste0(group, '\n(n = ', group_size, ')'))

# plot ...

p1 <- ggplot(data=filter(plot, study == 'ORC'))
p2 <- ggplot(data=filter(plot, study == 'SRC'))

s <- list(
  aes(x = rating, y = token, group = environment, col = environment, 
      fill = environment, linetype = environment),
  geom_hline(yintercept = 2.5),
  geom_vline(xintercept = 3.5),
  geom_smooth(method=lm, alpha = .2),
  theme_classic(),
  coord_cartesian(ylim = c(-.5, 5.5),
                  xlim = c(1, 6)),
  scale_x_continuous(name='Rating',
                     breaks = c(1, 3.5, 6)),
  scale_y_continuous(name='Tokens', 
                     breaks = c(0, 2.5, 5)),
  scale_linetype_manual(name = 'Environment',
                        values = c('dashed', 'solid', 'twodash'),
                        labels = c('Short', 'Long', 'Island')),
    scale_color_manual(name = 'Environment',
                       values = c('#619CFF', '#00BA38', '#F8766D'),
                       labels = c('Short', 'Long', 'Island')),
    scale_fill_manual(name = 'Environment',
                      values = c('#619CFF', '#00BA38', '#F8766D'),
                      labels = c('Short', 'Long', 'Island')),
  theme(text = element_text(size = 12),
        plot.margin = margin(2, 2, 2, 2),
        legend.margin = margin(2, 2, 2, 2),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "right"),
  facet_grid2(~group_label, axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ept_ajt_token_rating.png', width=6.5, height=2, dpi=600)

p2 + s
ggsave('plots/src_ept_ajt_token_rating.png', width=6.5, height=2, dpi=600)

# prep dataset ...

md1 <- plot %>%
  mutate(group = factor(group))

contrasts(md1$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md1$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md1$environment)
contrasts(md1$group)

md2 <- plot %>%
  filter(environment %in% c('Long', 'Island')) %>%
  mutate(environment = fct_drop(environment)) %>%
  mutate(group = factor(group))

summary(md2$environment)
summary(md2$group)

contrasts(md2$environment) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)
contrasts(md2$group) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)

contrasts(md2$environment)
contrasts(md2$group)

# fit negative binomial regression model ...

# ... orc model (final) ####

library(glmmTMB)
model <- glmmTMB::glmmTMB(token ~ environment * group * zscore + (1 | participant),
                  data = filter(md1, study == 'ORC'),
                  ziformula = ~0,
                  family = nbinom2)
summary(model)

check_convergence(model) # TRUE
check_singularity(model) # FALSE
check_overdispersion(model) # FALSE

# no warnings

# ... src model (final) ####

model <- glmmTMB::glmmTMB(token ~ environment * group * rating + (1 | participant),
                          data = filter(md2, study == 'SRC'),
                          ziformula = ~0,
                          family = nbinom2)
summary(model)

check_convergence(model) # TRUE
check_singularity(model) # FALSE
check_overdispersion(model) # FALSE

# no warnings

# ---------------------------------------------------------------------------- #
# plot: proficiency vs. consistent RP acceptance for L2ers ####
# ---------------------------------------------------------------------------- #

# plot by L2 proficiency

ajt_ind2 <- ajt_ind %>%
  filter(pairing %in% c('KE', 'ME')) %>%
  filter(is.na(proficiency) == FALSE) %>%
  mutate(proficiency = as.numeric(proficiency)) %>%
  mutate(acceptor = case_when(type == 'accept' ~ TRUE, TRUE ~ FALSE))

plot <- ajt_ind2 %>%
  group_by(study, pairing) %>%
  mutate(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in ', language, '\n(n = ', pairing_size, ')')) %>%
  group_by(study, pairing, pairing_label, environment, acceptor) %>%
  summarise(mean = mean(proficiency),
            lcl = my.lcl(proficiency),
            ucl = my.ucl(proficiency)) %>%
  ungroup()

p1 <- ggplot(data = filter(plot, study == 'ORC'))
p2 <- ggplot(data = filter(plot, study == 'SRC'))

s <- list(
  aes(x=environment, y=mean, group=acceptor, col=acceptor, shape=acceptor),
  geom_line(position = position_dodge(width = .4)),
  geom_point(size = 3, alpha = .65, position = position_dodge(width = .4)),
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2, linetype = 1, 
                position = position_dodge(width = .4)),
  theme_classic(),
  scale_x_discrete(name='Environment', 
                   limits = c('short', 'long', 'island'), 
                   labels = c('Short', 'Long', 'Island'),
  ),
  scale_y_continuous(name='Proficiency Score', limits=c(0, 50)),
  scale_colour_manual(name='Type', 
                      values=c('#40B0A6', '#E1BE6A'), 
                      limits = c('TRUE', 'FALSE'),
                      labels=c('RP Acceptor', 'RP Non-Acceptor')),
  scale_shape_manual(name='Type', 
                     values=c(15, 16), 
                     limits = c('TRUE', 'FALSE'),
                     labels=c('RP Acceptor', 'RP Non-Acceptor')),
  theme(
    text = element_text(size = 10),
    legend.position = 'right',
    plot.margin = margin(0, 0, 0, 1)),
  facet_manual(~pairing_label, design = 'AB', axes = 'all', remove_labels = FALSE)
)

p1 + s
ggsave('plots/orc_ajt_acceptor_proficiency.png', width=6.5, height=2.5, dpi=600)

p2 + s
ggsave('plots/src_ajt_acceptor_proficiency.png', width=6.5, height=2.5, dpi=600)

# ---------------------------------------------------------------------------- #
# model: proficiency vs. consistent RP acceptance for L2ers ####
# ---------------------------------------------------------------------------- #

# prep dataframe ...

md <- ajt_ind2 %>%
  select(study, pairing, participant, acceptor, environment, proficiency) %>%
  mutate_at(c('study', 'pairing', 'participant'), factor) %>%
  mutate(proficiency = scale(as.numeric(proficiency))) %>%
  mutate(proficiency = as.numeric(proficiency)) %>%
  mutate(environment = fct_relevel(environment, 'short', 'long', 'island'))

summary(md)

# apply deviation coding ...

contrasts(md$environment) <- contr.treatment(3) - matrix(rep(1/3, 6), ncol=2)
contrasts(md$pairing) <- contr.treatment(2) - matrix(rep(1/2, 2), ncol=1)

# view contrasts ...

contrasts(md$environment)
contrasts(md$pairing)

# ---------------------------------------------------------------------------- #
# ... orc model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model <- glmer(acceptor ~ environment * pairing * proficiency + 
                 (1 | participant), 
               data = filter(md, study == 'ORC'), family = binomial, 
               control = glmerControl(
                 optimizer = "bobyqa", 
                 optCtrl=list(maxfun=1e6)))
summary(model)
toc()
beep()

# check for warnings ...

performance::check_convergence(model) # TRUE
performance::check_singularity(model) # FALSE

# save model ...

model %>% write_rds('models/orc_ajt_acceptor_proficiency_model.rds')

# load model (if you have run and saved it previously) ...

model <- read_rds('models/orc_ajt_acceptor_proficiency_model.rds')

# model summary ...

summary(model)

# Parameter                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        -1.5321     0.4033  -3.799 0.000145 ***
# environment2                        0.4144     0.4008   1.034 0.301216    
# environment3                        0.8408     0.4117   2.042 0.041127 *  
# pairing2                            0.9668     0.6838   1.414 0.157403    
# proficiency                        -0.7641     0.3646  -2.096 0.036112 *  
# environment2:pairing2              -0.8930     0.8024  -1.113 0.265781    
# environment3:pairing2              -2.5940     0.8444  -3.072 0.002128 ** 
# environment2:proficiency            0.2154     0.3995   0.539 0.589719    
# environment3:proficiency            0.1531     0.4055   0.378 0.705693    
# pairing2:proficiency                0.3759     0.6850   0.549 0.583186    
# environment2:pairing2:proficiency  -0.8533     0.8021  -1.064 0.287369    
# environment3:pairing2:proficiency  -1.1432     0.8205  -1.393 0.163531 

# table ...

table <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  rename('SE' = 'Std. Error', 'z' = 'z value', 'p' = 'Pr(>|z|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'pairing2' = 'Pairing[ME]',
              'proficiency' = 'Proficiency'))

table %>% write_rds('tables/orc_ajt_acceptor_proficiency_model.rds')

# joint tests ...

joint_tests(model)

# model term                      df1 df2 F.ratio p.value
# environment                       2 Inf   2.117  0.1204
# pairing                           1 Inf   1.918  0.1660
# proficiency                       1 Inf   4.392  0.0361 *
# environment:pairing               2 Inf   4.871  0.0077 **
# environment:proficiency           2 Inf   0.150  0.8605
# pairing:proficiency               1 Inf   0.301  0.5832
# environment:pairing:proficiency   2 Inf   1.033  0.3559

# post-hoc tests ...

emtrends <- emtrends(model, ~ pairing * environment, 
                     var = 'proficiency', 
                     adjust = 'mvt')
test <- test(emtrends)
test

# pairing environment proficiency.trend    SE  df z.ratio p.value
# KE      short                  -1.408 0.708 Inf  -1.989  0.2139
# ME      short                  -0.366 0.522 Inf  -0.701  0.9645
# KE      long                   -0.766 0.643 Inf  -1.191  0.7252
# ME      long                   -0.578 0.527 Inf  -1.096  0.7884
# KE      island                 -0.683 0.613 Inf  -1.114  0.7765
# ME      island                 -0.785 0.541 Inf  -1.450  0.5383

# table ...

table <- test %>%
  rename('Pairing' = 'pairing', 'Environment' = 'environment', 
         'Proficiency Trend' = 'proficiency.trend', 'z' = 'z.ratio',
         'p' = 'p.value') %>%
  select(Pairing, Environment, everything()) %>%
  mutate_at(c('Pairing', 'Environment'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Pairing, Environment)

table %>% write_rds('tables/orc_ajt_acceptor_proficiency_posthoc.rds')

# plot ...

emmip(model, ~ environment | pairing, var = 'proficiency', cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

label <- md %>%
  mutate(group = case_when(pairing == 'KE' ~ 'KLE', 
                           pairing == 'ME' ~ 'MLE')) %>%
  group_by(study, group, pairing) %>%
  summarise(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in English\n(n = ', pairing_size, ')')) %>%
  filter(study == 'ORC') %>%
  select(pairing, pairing_label)

plot <- emmip(model, proficiency ~ environment | pairing, cov.reduce = range, CIs = TRUE,
              combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = FALSE) %>%
  left_join(label, by = 'pairing')

ggplot(plot, aes(x = proficiency, y = yvar, group = environment, col = environment, linetype = environment)) +
  theme_classic() +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = environment), alpha = .1, col = NA) +
  geom_line(lwd = .75) +
  scale_y_continuous(name='Logits') +
  scale_x_continuous(name='Proficiency (z-score)') +
  scale_linetype_manual(name = 'Environment',
                        values = c('longdash', 'solid', 'twodash'),
                        labels = c('Short', 'Long', 'Island')) +
  scale_color_manual(name = 'Environment',
                     values = c('#619CFF', '#00BA38', '#F8766D'),
                     labels = c('Short', 'Long', 'Island')) +
  scale_fill_manual(name = 'Environment',
                    values = c('#619CFF', '#00BA38', '#F8766D'),
                    labels = c('Short', 'Long', 'Island')) +
  theme(legend.key.size = unit(2,"line")) +
  facet_wrap(~pairing_label)

ggsave('plots/orc_ajt_rating_proficiency_posthoc.png', width=6.5, height=2.75, dpi=600)

# ---------------------------------------------------------------------------- #
# ... src model (final) ####
# ---------------------------------------------------------------------------- #

tic()
model <- glmer(acceptor ~ environment * pairing * proficiency + 
                 (1 | participant), 
               data = filter(md, study == 'SRC'), family = binomial, 
               control = glmerControl(
                 optimizer = "bobyqa", 
                 optCtrl=list(maxfun=1e6)))
summary(model)
toc()
beep()

# check for warnings ...

performance::check_convergence(model) # TRUE
performance::check_singularity(model) # FALSE

# save model ...

model %>% write_rds('models/src_ajt_acceptor_proficiency_model.rds')

# load model (if you have run and saved it previously) ...

model <- read_rds('models/src_ajt_acceptor_proficiency_model.rds')

# model summary ...

summary(model)

# Parameter                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                       -0.45303    0.19341  -2.342  0.01917 *  
# environment2                       1.02214    0.34227   2.986  0.00282 ** 
# environment3                       1.17872    0.34164   3.450  0.00056 ***
# pairing2                           1.63198    0.40252   4.054 5.03e-05 ***
# proficiency                       -0.38404    0.19284  -1.992  0.04642 *  
# environment2:pairing2              0.08279    0.66829   0.124  0.90141    
# environment3:pairing2             -1.09645    0.66687  -1.644  0.10014    
# environment2:proficiency           0.53491    0.34389   1.555  0.11983    
# environment3:proficiency           0.91127    0.34589   2.635  0.00842 ** 
# pairing2:proficiency               0.02627    0.37993   0.069  0.94486    
# environment2:pairing2:proficiency  0.30578    0.68375   0.447  0.65472    
# environment3:pairing2:proficiency -0.64394    0.67934  -0.948  0.34318  

# table ...

table <- summary(model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column('Parameter') %>%
  rename('SE' = 'Std. Error', 'z' = 'z value', 'p' = 'Pr(>|z|)') %>%
  mutate_if(is.character, str_replace_all, 
            c('environment2' = 'Environment[Long]',
              'environment3' = 'Environment[Island]',
              'pairing2' = 'Pairing[ME]',
              'proficiency' = 'Proficiency'))

table %>% write_rds('tables/src_ajt_acceptor_proficiency_model.rds')

# joint tests ...

joint_tests(model)

# model term                      df1 df2 F.ratio p.value
# environment                       2 Inf   7.291  0.0007 ***
# pairing                           1 Inf  16.550  <.0001 ***
# proficiency                       1 Inf   3.966  0.0464 *
# environment:pairing               2 Inf   2.405  0.0903 .
# environment:proficiency           2 Inf   3.471  0.0311 *
# pairing:proficiency               1 Inf   0.005  0.9449
# environment:pairing:proficiency   2 Inf   1.175  0.3089

# post-hoc tests ...

emtrends <- emtrends(model, ~ pairing * environment, 
                     var = 'proficiency', 
                     adjust = 'mvt')
test <- test(emtrends)
test

# pairing environment proficiency.trend    SE  df z.ratio p.value
# KE      short                  -0.936 0.465 Inf  -2.012  0.2309
# ME      short                  -0.797 0.364 Inf  -2.186  0.1568
# KE      long                   -0.554 0.400 Inf  -1.385  0.6481
# ME      long                   -0.109 0.345 Inf  -0.316  0.9997
# KE      island                  0.298 0.379 Inf   0.786  0.9623
# ME      island                 -0.207 0.338 Inf  -0.613  0.9892

# table ...

table <- test %>%
  rename('Pairing' = 'pairing', 'Environment' = 'environment', 
         'Proficiency Trend' = 'proficiency.trend', 'z' = 'z.ratio',
         'p' = 'p.value') %>%
  select(Pairing, Environment, everything()) %>%
  mutate_at(c('Pairing', 'Environment'), as.character) %>%
  mutate_if(is.character, str_replace_all, 
            c('short' = 'Short', 'long' = 'Long', 'island' = 'Island')) %>%
  select(-df) %>%
  mutate(Environment = fct_relevel(Environment, 'Short', 'Long', 'Island')) %>%
  arrange(Pairing, Environment)

table %>% write_rds('tables/src_ajt_acceptor_proficiency_posthoc.rds')

# plot ...

emmip(model, ~ environment | pairing, var = 'proficiency', cov.reduce = range, CIs = TRUE,
      combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = TRUE)

label <- md %>%
  mutate(group = case_when(pairing == 'KE' ~ 'KLE', 
                           pairing == 'ME' ~ 'MLE')) %>%
  group_by(study, group, pairing) %>%
  summarise(pairing_size = length(unique(participant))) %>%
  ungroup() %>%
  mutate(pairing_label = paste0(group, ' in English\n(n = ', pairing_size, ')')) %>%
  filter(study == 'ORC') %>%
  select(pairing, pairing_label)

plot <- emmip(model, proficiency ~ environment | pairing, cov.reduce = range, CIs = TRUE,
              combine = TRUE, comparisons = TRUE, adjust = 'mvt', plotit = FALSE) %>%
  left_join(label, by = 'pairing')

ggplot(plot, aes(x = proficiency, y = yvar, group = environment, col = environment, linetype = environment)) +
  theme_classic() +
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = environment), alpha = .1, col = NA) +
  geom_line(lwd = .75) +
  scale_y_continuous(name='Logits') +
  scale_x_continuous(name='Proficiency (z-score)') +
  scale_linetype_manual(name = 'Environment',
                        values = c('longdash', 'solid', 'twodash'),
                        labels = c('Short', 'Long', 'Island')) +
  scale_color_manual(name = 'Environment',
                     values = c('#619CFF', '#00BA38', '#F8766D'),
                     labels = c('Short', 'Long', 'Island')) +
  scale_fill_manual(name = 'Environment',
                    values = c('#619CFF', '#00BA38', '#F8766D'),
                    labels = c('Short', 'Long', 'Island')) +
  theme(legend.key.size = unit(2,"line")) +
  facet_wrap(~pairing_label)

ggsave('plots/src_ajt_rating_proficiency_posthoc.png', width=6.5, height=2.75, dpi=600)
