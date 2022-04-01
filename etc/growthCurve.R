library(tidyverse)
library(data.table)


Data <- read.csv('https://raw.githubusercontent.com/tlobnow/Resources/main/etc/growth%20curve.csv', sep = ';', na.strings = c('', ' ', 'NA'))

Data$WT <- Data$WT %>% str_replace_all(",", ".")
Data$WT.1 <- Data$WT.1 %>% str_replace_all(",", ".")
Data$WT.2 <- Data$WT.2 %>% str_replace_all(",", ".")
Data$WT.3 <- Data$WT.3 %>% str_replace_all(",", ".")
Data$WT.4 <- Data$WT.4 %>% str_replace_all(",", ".")
Data$S1.3 <- Data$S1.3 %>% str_replace_all(",", ".")
Data$S1.3.1 <- Data$S1.3.1 %>% str_replace_all(",", ".")
Data$S1.3.2 <- Data$S1.3.2 %>% str_replace_all(",", ".")
Data$S1.3.3 <- Data$S1.3.3 %>% str_replace_all(",", ".")
Data$S1.3.4 <- Data$S1.3.4 %>% str_replace_all(",", ".")

Data$WT <- as.numeric(Data$WT)
Data$WT.1 <- as.numeric(Data$WT)
Data$WT.2 <- as.numeric(Data$WT.2)
Data$WT.3 <- as.numeric(Data$WT.3)
Data$WT.4 <- as.numeric(Data$WT.4)
Data$S1.3 <- as.numeric(Data$S1.3)
Data$S1.3.1 <- as.numeric(Data$S1.3.1)
Data$S1.3.2 <- as.numeric(Data$S1.3.2)
Data$S1.3.3 <- as.numeric(Data$S1.3.3)
Data$S1.3.4 <- as.numeric(Data$S1.3.4)

Data_MT <- Data %>% dplyr::select(Day, WT) %>% 
  pivot_longer(names_to = "Temp", values_to = "WT", cols = c(WT)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)

Data_MT1 <- Data %>% dplyr::select(Day, WT.1) %>% 
  pivot_longer(names_to = "Temp", values_to = "WT.1", cols = c(WT.1)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)


Data_MT2 <- Data %>% dplyr::select(Day, WT.2) %>% 
  pivot_longer(names_to = "Temp", values_to = "WT.2", cols = c(WT.2)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)

Data_MT3 <- Data %>% dplyr::select(Day, WT.3) %>% 
  pivot_longer(names_to = "Temp", values_to = "WT.3", cols = c(WT.3)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)

Data_MT4 <- Data %>% dplyr::select(Day, WT.4) %>% 
  pivot_longer(names_to = "Temp", values_to = "WT.4", cols = c(WT.4)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)


Data_ST <- Data %>% dplyr::select(Day, S1.3) %>% 
  pivot_longer(names_to = "Temp", values_to = "S1.3", cols = c(S1.3)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)

Data_ST1 <- Data %>% dplyr::select(Day, S1.3.1) %>% 
  pivot_longer(names_to = "Temp", values_to = "S1.3.1", cols = c(S1.3.1)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)


Data_ST2 <- Data %>% dplyr::select(Day, S1.3.2) %>% 
  pivot_longer(names_to = "Temp", values_to = "S1.3.2", cols = c(S1.3.2)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)

Data_ST3 <- Data %>% dplyr::select(Day, S1.3.3) %>% 
  pivot_longer(names_to = "Temp", values_to = "S1.3.3", cols = c(S1.3.3)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)

Data_ST4 <- Data %>% dplyr::select(Day, S1.3.4) %>% 
  pivot_longer(names_to = "Temp", values_to = "S1.3.4", cols = c(S1.3.4)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  dplyr::ungroup() %>% 
  distinct(Day, .keep_all = T)


Data_WT <- full_join(Data_MT, Data_MT1) %>% filter(!is.na(WT) | !is.na(WT.1))
Data_WT <- full_join(Data_WT, Data_MT2) %>% filter(!is.na(WT) | !is.na(WT.1) | !is.na(WT.2))
Data_WT <- full_join(Data_WT, Data_MT3) %>% filter(!is.na(WT) | !is.na(WT.1) | !is.na(WT.2) | !is.na(WT.3))
Data_WT <- full_join(Data_WT, Data_MT4) %>% filter(!is.na(WT) | !is.na(WT.1) | !is.na(WT.2) | !is.na(WT.3) | !is.na(WT.4))

Data_S1.3 <- full_join(Data_ST, Data_ST1) %>% filter(!is.na(S1.3) | !is.na(S1.3.1))
Data_S1.3 <- full_join(Data_S1.3, Data_ST2) %>% filter(!is.na(S1.3) | !is.na(S1.3.1) | !is.na(S1.3.2))
Data_S1.3 <- full_join(Data_S1.3, Data_ST3) %>% filter(!is.na(S1.3) | !is.na(S1.3.1) | !is.na(S1.3.2) | !is.na(S1.3.3))
Data_S1.3 <- full_join(Data_S1.3, Data_ST4) %>% filter(!is.na(S1.3) | !is.na(S1.3.1) | !is.na(S1.3.2) | !is.na(S1.3.3) | !is.na(S1.3.4))

Data <- full_join(Data_WT, Data_S1.3)

Data_pivot <- Data %>% dplyr::select(Day, WT.1, WT.2, WT.3, WT.4, S1.3, S1.3.1, S1.3.2, S1.3.3, S1.3.4) %>% 
  pivot_longer(names_to = "Mouse", values_to = "Parasitemia", cols = c(WT.1, WT.2, WT.3, WT.4, S1.3, S1.3.1, S1.3.2, S1.3.3, S1.3.4)) %>% 
  dplyr::arrange(Day) %>% 
  dplyr::group_by(Day) %>%
  dplyr::filter(!is.na(Parasitemia))

Data_vis <- Data_pivot



Data_vis <- Data_vis %>% mutate(Type = case_when(Mouse == 'WT' ~ 'WT',
                                               Mouse == 'WT.1' ~ 'WT',
                                               Mouse == 'WT.2' ~ 'WT',
                                               Mouse == 'WT.3' ~ 'WT',
                                               Mouse == 'WT.4' ~ 'WT',
                                               Mouse == 'S1.3' ~ 'sera1-3',
                                               Mouse == 'S1.3.1' ~ 'sera1-3',
                                               Mouse == 'S1.3.2' ~ 'sera1-3',
                                               Mouse == 'S1.3.3' ~ 'sera1-3',
                                               Mouse == 'S1.3.4' ~ 'sera1-3'))

RT_sum <- Data_vis %>% group_by(Day, Type) %>% filter(Mouse != 'S1.3.4') %>% summarise(Mean = mean(Parasitemia, na.rm=TRUE),
                                                         SD = sd(Parasitemia, na.rm = TRUE))
Data_vis <- full_join(Data_vis, RT_sum)

Data_vis %>%
  ggplot(aes(x = Day, col = Type)) +
  #geom_smooth(aes(y = Mean)) +
  #geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, col = Type), width=.1) +
  geom_ribbon(aes(ymin=Mean-SD, ymax=Mean+SD, fill = Type), col = 'white',  alpha = 0.4) +
  geom_line(aes(y = Mean)) +
  theme_classic()
  




