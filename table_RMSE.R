# Generate RMSE tables 
load("rmse.RData")

row.idx = c(3, 4, 1, 2)

table.data = cbind(rmse_100[row.idx, ], rmse_500[row.idx, ], rmse_1000[row.idx, ])

# rescale by 100, round by 2 digits

table1 = round(table.data * 100, 2)
library(xtable)
xtable(table1)

# Generate number of multiple \hat{t}'s 

library(readxl)
table_multiple <- read_excel("table_multiple.xlsx", 
                               +     col_names = FALSE)
xtable(table_multiple)
