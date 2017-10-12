#ASSOCIATION RULE MINING OF STB PATIENT META DATA AND SEQUENCING DATA IN R
#USING THE ARULES PACKAGE

#FILE CONTAINS METADATA & SEQUENCING DATA FOR 25 STB STRAINS
#HEADINGS INCLUDE:
#1. ID
#2. Genome length in Mb
#3. Number of CDS
#4. Number of SNPs
#5. Year of Isolation
#6. Patient Nationality
#7. Patient Gender
#8. Site of infection
#9. Strain type
#10. Patient Age

#NUMERICAL VALUES WERE CONVERTED TO CATEGORICAL RANGES
#AIM TO ASSOCIATE VARIABLES TO PIN POINT CORRELATIONS BETWEEN PATIENTS 
#AND STRAIN TYPE TO DEDUCE ANY STB INFECTION PROFILES
##################################################################

#1. LOAD ARULES AND INPUT METADATA FILE - NEED TO HIGHLIGHT N/A INFO

required_data <- read.table("25_metadata4arules.txt",sep="\t", na.strings="N/A",h=T,stringsAsFactors=TRUE)
required_data <- as(required_data,"transactions")

##################################################################

#2. CREATE CLUSTER AND FREQUENCY PLOTS TO GET A VIEW OF THE DATA

#CLUSTER PLOT
d <- dissimilarity(sample(required_data, 25), method = "phi", which = "items")
d[is.na(d)] <- 1
plot(hclust(d), cex=.5)

#MOST FREQUENT FACTORS < 10%
itemFrequencyPlot(required_data,support=0.1,cex.names=0.8)

##################################################################

#3A. RUN ARULES WITH A MINIMUM LENGTH OF 5, #SUPPORT OF 0.15 AND A CONFIDENCE OF 0.6
#SET RHS TO STRAIN TYPE

rules <- apriori(required_data, parameter=list(minlen=5,supp=0.15
,conf=0.6), appearance=list(rhs=c("STRAIN=STB-A","STRAIN=STB-D","STRAIN=STB-L","STRAIN=STB-F","STRAIN=STB-H","STRAIN=STB-E","STRAIN=STB-J","STRAIN=STB-G","STRAIN=STB-I","STRAIN=STB-N","STRAIN=STB-K"), default="lhs"))

#3B. RUN ARULES WITH A MINIMUM LENGTH OF 5, #SUPPORT OF 0.15 AND A CONFIDENCE OF 0.6
#SET RHS TO NATIONALITY

rules <- apriori(required_data, parameter=list(minlen=5,supp=0.15,conf=0.6), appearance=list(rhs=c("NATIONALITY=FR","NATIONALITY=DJ","NATIONALITY=ETH"), default="lhs"))

##################################################################

#4A. SORT RESULTS TABLE BY "LIFT" - HIGHER THE VALUE
#THE STRONGER THE ASSOCIATION 

rules.sorted <- sort(rules, by="lift")

#4B. HAVE A QUICK LOOK AT THE TOP 6 RULES

inspect(head(rules.sorted))

##################################################################

#A5. PPLY A FISHER'S EXTACT TEST TO VALIDATE FINDINGS

interestMeasure(rules.sorted, measure="fishersExactTest", transactions=required_data)
