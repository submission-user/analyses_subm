library(psychotree)
?SPISA
#A subsample from the general knowledge quiz “Studentenpisa” conducted online by the German 
# weekly news magazine SPIEGEL. The data contain the quiz results from 45 questions as well 
# as sociodemographic data for 1075 university students from Bavaria.
# spisa
# matrix with 0/1 results from 45 questions in the quiz (indicating wrong/correct answers).
# 
# gender
# factor indicating gender.
# 
# age
# age in years.
# 
# semester
# numeric indicating semester of university enrollment.
# 
# elite
# factor indicating whether the university the student is enrolled in has been granted “elite” 
# status by the German “excellence initiative”.
# 
# spon
# ordered factor indicating frequency of accessing the SPIEGEL online (SPON) magazine.
#
#
# Details
# An online quiz for testing one's general knowledge was conducted by the German weekly news 
# magazine SPIEGEL in 2009. Overall, about 700,000 participants answered the quiz and a set of 
# sociodemographic questions. The general knowledge quiz consisted of a total of 45 items from 
# five different topics: politics, history, economy, culture and natural sciences. For each topic,
# four different sets of nine items were available, that were randomly assigned to the participants. 
# A thorough analysis and discussion of the original data set is provided in Trepte and Verbeet (2010).
# 
# Here, we provide the subsample of university students enrolled in the federal state of Bavaria, 
# who had been assigned questionnaire number 20 (so that all subjects have answered the same set of 
#                                                items). Excluding all incomplete records, this 
# subsample contains 1075 observations.
# 
# The data are analyzed in Strobl et al. (2010), whose analysis is replicated in vignette("raschtree", package = "psychotree").
# 
# The full list of items in questionnaire 20 is given below.
# 
# Politics:
# Who determines the rules of action in German politics according to the constitution? – The Bundeskanzler (federal chancellor).
# What is the function of the second vote in the elections to the German Bundestag (federal parliament)? – It determines the allocation of seats in the Bundestag.
# How many people were killed by the RAF (Red Army Faction)? – 33.
# Where is Hessen (i.e., the German federal country Hesse) located? – (Indicate location on a map.)
# What is the capital of Rheinland-Pfalz (i.e., the German federal country Rhineland-Palatinate)? – Mainz.
# Who is this? – (Picture of Horst Seehofer.)
# Which EU institution is elected in 2009 by the citizens of EU member countries? – European Parliament.
# How many votes does China have in the UNO general assembly? – 1.
# Where is Somalia located? – (Indicate location on a map.)
# 
# History:
# The Roman naval supremacy was established through... – ... the abolition of Carthage.
# In which century did the Thirty Years' War take place? – The 17th century.
# Which form of government is associated with the French King Louis XIV? – Absolutism.
# What island did Napoleon die on in exile? – St. Helena.
# How many percent of the votes did the NSDAP receive in the 1928 elections of the German Reichstag? – About 3 percent.
# How many Jews were killed by the Nazis during the Holocaust? – About 6 Million.
# Who is this? – (Picture of Johannes Rau, former German federal president.)
# Which of the following countries is not a member of the EU? – Croatia.
# How did Mao Zedong expand his power in China? – The Long March.
# 
# Economy:
#   Who is this? – (Picture of Dieter Zetsche, CEO of Mercedes-Benz.)
# What is the current full Hartz IV standard rate (part of the social welfare) for adults? – 351 Euro.
# What was the average per capita gross national product in Germany in 2007? – About 29,400 Euro.\ What is a CEO? – A Chief Executive Officer.
# What is the meaning of the hexagonal “organic” logo? – Synthetic pesticides are prohibited.
# Which company does this logo represent? – Deutsche Bank.
# Which German company took over the British automobile manufacturers Rolls-Royce? – BMW.
# Which internet company took over the media group Time Warner? – AOL.
# What is the historic meaning of manufacturies? – Manufacturies were the precursors of industrial mass production.
# Culture:
#   Which painter created this painting? – Andy Warhol.
# What do these four buildings have in common? – All four were designed by the same architects.
# Roman numbers: What is the meaning of CLVI? – 156.
# What was the German movie with the most viewers since 1990? – Der Schuh des Manitu.
# In which TV series was the US president portrayed by an African American actor for a long time? – 24.
# What is the name of the bestselling novel by Daniel Kehlmann? – Die Vermessung der Welt (Measuring The World).
# Which city is the setting for the novel ‘Buddenbrooks’? – Lübeck.
# In which city is this building located? – Paris.
# Which one of the following operas is not by Mozart? – Aida.
# 
# Natural sciences:
#   Why does an ice floe not sink in the water? – Due to the lower density of ice.
# What is ultrasound not used for? – Radio.
# Which sensory cells in the human eye make color vision possible? – Cones.
# What is also termed Trisomy 21? – Down syndrome.
# Which element is the most common in the Earth's atmosphere? – Nitrogen.
# Which kind of tree does this leaf belong to? – Maple.
# Which kind of bird is this? – Blackbird.
# Where is the stomach located? – (Indicate location on a map of the body.)
# What is the sum of interior angles in a triangle? – 180 degrees.
# 
# References
# Strobl C, Kopf J, Zeileis A (2015). Rasch Trees: A New Method for Detecting Differential 
# Item Functioning in the Rasch Model. Psychometrika, 80(2), 289–316. doi: 10.1007/s11336-013-9388-3
# 
# SPIEGEL Online (2009). Studentenpisa – Alle Fragen, alle Antworten. In German. Accessed 2010-10-26. 
# https://www.spiegel.de/lebenundlernen/uni/studentenpisa-alle-fragen-alle-antworten-a-620101.html
# 
# Trepte S, Verbeet M (2010). Allgemeinbildung in Deutschland – Erkenntnisse aus dem 
# SPIEGEL-Studentenpisa-Test. ISBN 978-3-531-17218-7. VS Verlag, Wiesbaden.




library(sirt)
# n=565, 11 0/1-items, 2 0/1-covs, 1 cont cov
?data.pisaMath

# Description
# This is an example PISA dataset of mathematics items. The dataset contains 565 students on 11 items.
# 
# Usage
# data(data.pisaMath)
# Format
# The dataset is a list. The list element data contains the dataset with the demographical variables 
# student ID (idstud), school ID (idschool), a dummy variable for female students (female), socioeconomic 
# status (hisei) and migration background (migra). The remaining variables (starting with M in the name) 
# are the mathematics items. 
# The item metadata are included in the list element item which contains item name (item) and the testlet 
# label (testlet). An item not included in a testlet is indicated by NA.


# n=623, 12 0/1-items, 2 0/1-covs, 1 cont cov
?data.pisaRead

# Description
# This is an example PISA dataset of reading items. The dataset contains 623 students on 12 items.
# 
# Usage
# data(data.pisaRead)
# Format
# The dataset is a list. The list element data contains the dataset with the demographical 
# variables student ID (idstud), school ID (idschool), a dummy variable for female students 
# (female), socioeconomic status (hisei) and migration background (migra). The remaining variables
# (starting with R in the name) are the reading items. 
# The item metadata are included in the list element item which contains item name (item), testlet 
# label (testlet), item format (ItemFormat), text type (TextType) and text aspect (Aspect).


# n=345, 25 0/1-items, 1 0/1-covs, 1 cont cov
?data.timss

# Description
# This datasets contains TIMSS mathematics data from 345 students on 25 items.
# 
# Usage
# data(data.timss)
# Format
# This dataset is a list. data is the dataset containing student ID (idstud), 
# a dummy variable for female (girl) and student age (age). The following variables 
 #(starting with M in the variable name are items.
