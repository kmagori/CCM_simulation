#this is the R file that weaves in-line R code chunks found in the input ODF into text and graphics of the output ODF	
library(odfWeave)

odfWeave(file="CCM_simulation.odt",dest="test.odt",workDir="temp",control=odfWeaveControl(zipCmd=c("7z a $$file$$","7z x $$file$$")))