# Setup & load data

library(tidyverse)
setwd('/home/marten/AI7/MasterProject/extra/behav/')

orig <- read_csv('behavior_data.csv') %>%
  select(-`...1`, -fully_correct, -semi_correct, -any_correct) %>%
  filter(stim_type == 'color')

acc <- orig %>%
  group_by(session) %>%
  summarise(accuracy=mean(correct)) %>%
  arrange(desc(accuracy))

dat <- orig %>%
  pivot_longer(starts_with('image_num_'), names_pattern='image_num_p([12])_resp_1', names_to='subj', values_to='image_num') %>%
  inner_join(read_tsv('colorsUsed_2019-09-30_193104.dat'), by=c('image_set'='trial_set', 'image_num'='image')) %>%
  pivot_longer(c('unfixed_color', 'fixed_color1', 'fixed_color2'), names_to='color_type', values_to='color')


# Plot sessions

plot_session <- function (sess) {
  curdat <- dat %>%
    filter(session == sess) %>%
    group_by(subj, color_type, color) %>%
    summarise(count=n())

  ggplot(curdat, aes(color, count, fill=color, alpha=color_type)) +
    geom_bar(stat='identity', color='black') +
    facet_wrap(vars(paste('subj', subj))) +
    coord_flip() +
    scale_fill_manual(values=unique(curdat$color), guide='none') +
    scale_alpha_discrete(range=c(1/3, 1), guide=guide_legend(reverse=TRUE)) +
    theme(legend.position='bottom')
}

for (i in 1:nrow(acc)) {
  row <- acc[i,]
  plot_session(row$session) +
    labs(title=paste('Occurence of colors in chosen images: session', row$session, '; accuracy', round(row$accuracy, 2)))
  ggsave(str_c('out/color_freq_', i, '.png'))
}


# Plot average

curdat <- dat %>%
  group_by(color_type, color) %>%
  summarise(count=n())
ggplot(curdat, aes(color, count, fill=color, alpha=color_type)) +
  geom_bar(stat='identity', color='black') +
  coord_flip() +
  scale_fill_manual(values=unique(curdat$color), guide='none') +
  scale_alpha_discrete(range=c(1/3, 1), guide=guide_legend(reverse=TRUE)) +
  theme(legend.position='bottom')

dat %>%
  group_by(color_type, color, session) %>%
  summarise(count=n()) %>%
  mutate(random=runif(n())) %>%
  arrange(-count, random) %>%
  group_by(color_type, session) %>%
  filter(row_number() == 1) %>%
  ggplot(aes(color_type, fill=color)) +
    geom_bar() +
    coord_flip() +
    facet_wrap(vars(sprintf('session %02d', session))) +
    scale_fill_manual(values=unique(curdat$color), guide='none') +
    scale_x_discrete(limits=rev) +
    # thanks https://stackoverflow.com/a/1331400 !
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(str_c('out/summary.pdf'), width=18, height=12, units='cm')
