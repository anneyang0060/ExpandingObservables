setwd('C:/Users/annie/Downloads/accepted_2007_to_2018Q4.csv')
 
### exclude unrelevant varaibles
data<- read.csv('accepted_2007_to_2018Q4.csv')

# id
id_del <- c(1651666,1654417)
data <- data[-id_del, ]
id_num <- as.numeric(data$id)
id_na <- which(is.na(id_num))
for(i in id_na){
  n <- nchar(data$id[i])
  # print(i)
  id_num[i] <- as.numeric(substr(data$id[i], 39, n))
}
data$id <- id_num

# exclude unrelevant variables
excluded_vars <- c(
  # ---- Identifiers and metadata ----
  "member_id", "url", "desc", "title",

  # ---- Loan outcome or post-origination variables (leakage) ----
  "loan_status", "out_prncp", "out_prncp_inv",
  "total_pymnt", "total_pymnt_inv", "total_rec_prncp",
  "total_rec_int", "total_rec_late_fee", "recoveries",
  "collection_recovery_fee", "last_pymnt_d", "last_pymnt_amnt",
  "next_pymnt_d", "last_credit_pull_d",

  # ---- Hardship and settlement (post-loan events) ----
  "hardship_flag", "hardship_type", "hardship_reason",
  "hardship_status", "deferral_term", "hardship_amount",
  "hardship_start_date", "hardship_end_date",
  "payment_plan_start_date", "hardship_length",
  "hardship_dpd", "hardship_loan_status",
  "orig_projected_additional_accrued_interest",
  "hardship_payoff_balance_amount",
  "hardship_last_payment_amount",
  "debt_settlement_flag", "debt_settlement_flag_date",
  "settlement_status", "settlement_date", "settlement_amount",
  "settlement_percentage", "settlement_term",#
  # ---- Policy or constant fields ----
  "policy_code", "zip_code" 
)

excluded_idx <- c()
for(v in excluded_vars){
  j <- which(names(data)==v)
  excluded_idx <- c(excluded_idx,j)
}
data <- data[,-excluded_idx]

# order
order_id <- order(data$id)
data <- data[order_id, ]

#### Ordinal Encoding
# term
data$term[which(data$term==" 36 months")] <- 3
data$term[which(data$term==" 60 months")] <- 5
data$term <- as.numeric(data$term)

# grade
grade_num <- as.numeric(factor(data$grade, levels = c("A", "B", "C", "D", "E", "F", "G", "") ))
data$grade <- grade_num
lev <- sort(unique(data$sub_grade))
lev <- lev[-1]
lev <- c(lev, "")
grade_num <- as.numeric(factor(data$sub_grade, levels = lev ))
data$sub_grade <- grade_num
rm(grade_num)

# emp_length
lev <- sort(unique(data$emp_length))
num <- as.numeric(factor(data$emp_length, levels = lev ))
data$emp_length <- num

# pymnt_plan
lev <- sort(unique(data$pymnt_plan))
num <- as.numeric(factor(data$pymnt_plan, levels = lev ))
num[which(num<=2)] <- 0
num[which(num==3)] <- 1
data$pymnt_plan <- num

# verification status
num <- as.numeric(factor(
  data$verification_status, 
  levels = c("", "Not Verified", "Source Verified", "Verified"), 
  ordered = TRUE
))
data$verification_status <- num
num <- as.numeric(factor(
  data$verification_status_joint, 
  levels = c("", "Not Verified", "Source Verified", "Verified"), 
  ordered = TRUE
))
data$verification_status_joint <- num



# date
library(lubridate)

date1 <- data$issue_d
date0 <- data$earliest_cr_line
date1 <- decimal_date(my(date1))
date0 <- decimal_date(my(date0))
id_del <- which(is.na(date1-date0))
data$credit_history_length <- date1-date0
data <- data[-id_del,]
rm(date0,date1)
v <- which(names(data)=='earliest_cr_line')
data <- data[, -v]

dates <- data$issue_d
dates <- decimal_date(my(dates))
dates <- dates - min(dates)
data$issue_d <- dates

dates <- data$sec_app_earliest_cr_line
dates <- decimal_date(my(dates))-2007
data$sec_app_earliest_cr_line <- dates

### binary
data$initial_list_status <- ifelse(data$initial_list_status == "W", 1, 0)
data$application_type <- ifelse(data$application_type == "JOINT", 1, 0)
data$disbursement_method <- ifelse(data$disbursement_method == "DirectPay", 1, 0)

### categorical
top_jobs <- names(sort(table(data$emp_title), decreasing = TRUE)[1:51])[-1]
data$emp_title_simplified <- ifelse(data$emp_title %in% top_jobs, data$emp_title, "Other")
v <- which(names(data)=='emp_title')
data <- data[,-v]

manager = c("Manager",  "Project Manager", "Office Manager", "General Manager",
            "Operations Manager","Sales Manager","manager","Assistant Manager",
            "Store Manager", "Program Manager","Branch Manager", "Account Manager")
supervisor = c( "supervisor", "Supervisor ")
assistant = c("Administrative Assistant", "Executive Assistant")
owner = c("Owner", "owner")
president = c("Vice President","President")
teacher = c("Teacher","teacher")
nurse = c("Nurse","Registered Nurse","RN", "Registered nurse")
driver = c("Truck Driver","Driver","driver", "truck driver")
sales = c("Sales", "sales")
engineer = c("Engineer",  "Software Engineer")

data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% manager, "manager", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% supervisor, "supervisor", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% assistant, "assistant", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% owner, "owner", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% president, "president", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% nurse, "nurse", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% driver, "driver", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% sales, "sales", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% teacher, "teacher", data$emp_title_simplified)
data$emp_title_simplified <- ifelse(data$emp_title_simplified %in% sales, "sales", data$emp_title_simplified)


data$r1 <- data$funded_amnt/data$loan_amnt
data$r2 <- data$funded_amnt_inv/data$loan_amnt
data <- subset(data, select = -c(funded_amnt, funded_amnt_inv))

p <- ncol(data)
for(j in 2:p){
  if(typeof(data[,j])=="character") next
  idx <- which(!is.na(data[,j]))
  x <- data[idx,j]
  gap <- max(x)-min(x)
  if(gap>3){data[,j] <- data[,j]/gap}
}

# data <- read.csv('P2P_full.csv')
cat_vars <- c("home_ownership", "purpose", "addr_state", 'emp_title_simplified')
data <- dummy_cols(data, select_columns = cat_vars, 
                       remove_selected_columns = TRUE, remove_first_dummy = TRUE)

write.csv(data, 'P2P_full.csv')

### split to sub-datasets
K <- 2261
n <- c(rep(1000,K-1),639)
N_k <- 0
for(k in 1:K){
  idx <- (N_k + 1):(N_k + n[k])
  data_cur <- data[idx, ]
  file_cur <- paste0('P2Pdata_',k,'.csv')
  write.csv(data_cur, file_cur)
  N_k <- N_k + n[k]
}

pk <- list()
for(k in 1:K){
  print(k)
  file_cur <- paste0('P2Pdata_',k,'.csv')
  data_cur <- read.csv(file_cur)
  n <- nrow(data_cur)
  p <- ncol(data_cur)
  n_j <- sapply(1:p, function(j){sum(!is.na(data_cur[,j]))})
  p_obs <- which(n_j >= n*0.9)
  pk <- c(pk, list(p_obs))
}

p_mono <- list()
for(k in 1:K){
  p_obs <-1:p
  for (i in k:K) {
    print(paste0('k=',k,': ',i))
    p_obs <- intersect(p_obs, pk[[i]])
  }
  p_mono <- c(p_mono, list(p_obs))
}

for(k in 2:72){
  p_mono[[k]] <- p_mono[[1]]
}


for(k in 1:K){
  print(k)
  file_cur <- paste0('P2Pdata_',k,'.csv')
  data_cur <- read.csv(file_cur)
  data_cur <- data_cur[, p_mono[[k]]]
  data_cur <- data_cur[, -c(1,2)]
  write.csv(data_cur, file_cur)
}

L <- c()
for (k in 1:K) {
  L <- c(L, length(p_mono[[k]]))
}
idx <- 1:3
L_q <- unique(L)
for(ll in 4:length(L_q)){
  idx <- c(idx, min(which(L==L_q[ll])))
}
feature_set0 <- c()
for(k in idx){
  file_cur <- paste0('P2Pdata_',k,'.csv')
  data_cur <- read.csv(file_cur)
  feature_set1 <- names(data_cur)
  new_feature <- setdiff(feature_set1,feature_set0)
  feature_set0 <- feature_set1
  print(new_feature)
}
rm(list=ls())
gc()




