#####
#
# Example data (sirt pisaMath used in submission)
#
#####

# install.packages("difR")
# install.packages("difNLR")
# install.packages("ShinyItemAnalysis")
# install.packages("TAM")
# install.packages("sirt")
# install.packages("lordif")

load("datasets.RData")

# data from sirt


library(sirt)
# n=565, 11 0/1-items, 2 0/1-covs, 1 cont cov
# very well suited
data(data.pisaMath)
?data.pisaMath
names(data.pisaMath$data)

# n=623, 12 0/1-items, 2 0/1-covs, 1 cont cov
# very well suited
data(data.pisaRead)
?data.pisaRead
View(data.pisaRead)

# n=345, 25 0/1-items, 1 0/1-covs, 1 cont cov
# very well suited
data(data.timss)
?data.timss
View(data.timss)




# data from DIFlasso

library(psychotree)
data(SPISA)
?SPISA
names(SPISA)


# data from difR

library(difR)
data(verbal)
?verbal
names(verbal)
# The data set is called “verbal” and consists of 26 columns.
# The first 24 columns refer to the items, the 25th
# column (labeled “Anger”) corresponds to the trait anger
# score (Spielberger, 1988) of each subject, and the 26th
# column (labeled “Gender”) contains the group membership,
# with 0 and 1 entries for female and male respondents,
# respectively.


# data from difNLR
# n= 1407, 484 male0, 20 0/1-items, only one covariate
library(difNLR)
data(MSATB)
?MSATB
names(MSATB)

# data from ShinyItemAnalysis
# n=669, 20 items 0/1, two covariates 0/1

library(ShinyItemAnalysis)
data(HCI)
?HCI
names(HCI)

 #The first data set was collected during the final validation
# of the HCI (McFarland et al., 2017) and illustrates that finding
# a difference in total score between two groups does not 
# necessarily indicate item bias. The HCI is a 20-item 
# multiple-choice instrument designed to measure undergraduate 
# student understanding of homeostasis in physiology. The HCI was 
# validated with a sample of 669 undergraduate students, out of 
# whom 246 identified themselves as men, 405 identified themselves 
# as women, and the rest did not respond to this question 
# (McFarland et al., 2017). While the overall sample of 669 students
# is large, we knew that the sample sizes of the two subgroups might
# be small enough (ns < 500; in each group n < 500) that IRT models
# would be underpowered.

# data from TAM

library(TAM)

# 11 items scored, 3 covariates
# exclude items 3,4,7,8 , max score 2
data(data.timssAusTwn.scored, package = "TAM")
?data.timssAusTwn
names(data.timssAusTwn)

# n=1500, 157 items, three 0/1-covariates
# exlude non 0/1-items
# find description in Wu, Adams 2007
data(data.cqc05)
?data.cqc05
names(data.cqc05)

# n= 6371, 14 0/1-items, two 0/1-covariates (coded 1/2)
data(data.fims.Aus.Jpn.scored)
?data.fims.Aus.Jpn.scored
names(data.fims.Aus.Jpn.scored)


# data from lordif five category Likert

library(lordif)
data(Anxiety)
?Anxiety
names(Anxiety)



