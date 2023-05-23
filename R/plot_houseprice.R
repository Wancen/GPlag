# Load necessary libraries
library(ggplot2)
library(tidyverse)
baseline<-read.csv(file = "data/baselineprice2.csv")

baseline$county <- sub(",.*", "", baseline$region_name) # extract text before the comma
baseline$state <- sub(".*\\b([A-Z]{2})\\b.*", "\\1", baseline$region_name) # extract two capital letters surrounded by word boundaries

baseline %>% filter(state =="NC") -> res.price
res.price$county[which(res.price$county == "Durham")] = "Durham-Chapel Hill"
res.price$county[which(res.price$county == "Raleigh")] = "Raleigh-Cary"
res.price$county[which(res.price$county == "Greensboro")] = "Greensboro-High Point"
res.price$county[which(res.price$county == "Hickory")] = "Hickory-Lenoir-Morganton"


income<-read.csv(file = "data/personal_income.csv")%>% `colnames<-`(c("GeoFips", "GeoName","value"))
income$county <- sub(",.*", "", income$GeoName) # extract text before the comma
income$state <- sub(".*\\b([A-Z]{2})\\b.*", "\\1", income$GeoName) # extract two
income_nc <- income %>% filter(state =="NC")
income_price <- res.price %>% inner_join(income_nc)
corr_gplag <- cor(-income_price$alpha2, income_price$value, method = "spearman")
corr_tlcc <- cor(income_price$ccfcoef, income_price$value, method = "spearman")
corr_dtw <- cor(-income_price$dtwloss, income_price$value, method = "spearman")
corr_softdtw <- cor(-income_price$softdtwloss, income_price$value, method = "spearman")
corr_softdtwdiv <- cor(-income_price$softdtwdiv, income_price$value, method = "spearman")

corr_gplag <- cor(-income_price$alpha2, income_price$value, method = "kendall")
corr_tlcc <- cor(income_price$ccfcoef, income_price$value, method = "kendall")
corr_dtw <- cor(-income_price$dtwloss, income_price$value, method = "kendall")
corr_softdtw <- cor(income_price$softdtwloss, income_price$value, method = "kendall")
corr_softdtwdiv <- cor(-income_price$softdtwdiv, income_price$value, method = "kendall")

income_price$county <- factor(income_price$county, levels = income_price$county[order(income_price$value, decreasing = TRUE)])
income_price <- income_price %>%
  mutate(value_rank = rank(desc(value)),
         alpha_mle_rank = rank(alpha2),
         tlcc_rank = rank(desc(ccfcoef)),
         dtw_rank = rank(desc(dtwloss)),
         softdtw_rank = rank(softdtwloss),
         softdtwdiv_rank = rank(softdtwdiv))

income_rank <- income_price %>% select(county, value_rank:softdtwdiv_rank) %>% pivot_longer(cols = value_rank:softdtwdiv_rank,names_to = "method")
income_rank$method <- factor(income_rank$method, levels = c("value_rank"    , "alpha_mle_rank", "tlcc_rank"  ,    "dtw_rank","softdtw_rank", "softdtwdiv_rank"))

phouse <- ggplot(income_rank)+
  geom_linerange(aes(x = county, ymin = 0, ymax = value, colour = method),
                 position = position_dodge(width = 0.5),  linewidth = 1)+
  geom_point(aes(x = county, y = value, colour = method),
             position = position_dodge(width = 0.5),size = 2.5)+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.8),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("Rank") +
  scale_color_manual(values = c("red", "blue","green","grey","pink", "yellow"),
                     labels = c("Personal Income","GPlag","TLCC", "DTW","soft-DTW","soft-DTW divergence")) +
  guides(colour = guide_legend(title = "Methods"))
jpeg(file="plots/houseprice.jpg",width = 15, height = 4,units = "in",res=450)
phouse +
  theme(text = element_text(size = 14),
        plot.tag = element_text(size = 16)
  )
dev.off()



