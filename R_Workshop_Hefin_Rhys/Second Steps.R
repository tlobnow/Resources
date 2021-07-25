library(ggplot2)
library(dplyr)
library(mosaic)
library(car)
library(phia)


# reading in data ----
 pokemon <- read.csv("pokemon.csv")
    pokemon
  dim(pokemon)  
  head(pokemon)
  tail(pokemon)
  str(pokemon)  
  summary(pokemon)
  
# plotting ----
  # under name of variable --> Variable on X-axis
  # next to name of variable --> Variable on Y-Axis
    plot(pokemon)  
    plot(pokemon[ ,3:10])
  pokemon$Type.I <- as.factor(pokemon$Type.I) #because R changed stuff and now it would stop at zeros...
    plot(pokemon[ ,"Type.I"], pokemon[, "Atk"])


# t test ----  
  ?t.test
  
  psychic <-  pokemon[pokemon$Type.I == "Psychic", "Atk"]
              psychic  
  rock <-     pokemon[pokemon$Type.I == "Rock", "Atk"]
              rock    
  
  t.test(psychic, rock) #R doesn't assume that you have variables that are equal per default (FALSE)
  t.test(psychic, rock, var.equal = T)  # if you want to specifically request equal variables --> set this TRUE
                                        # leads to slightly different p-values
  # one-tailed t test
    t.test(psychic, rock, alternative = "less") # by default alternative is two-sided, but you can turn it one-tailed,
                                        # just use the side you want with "less than" or "greater than" other group, dep on order you put them
                                        # gives exactly half of previous p-value (--> one sided)
  # paired t test
    t.test(psychic[1:13], rock, paired = T) # must skip last two of psychic, since there is nothing to compare to in rock group (simply less)
     length(psychic)
      length(rock)  
  
  # for non-normal data --> Mann-Whitney u-test (in R weirdly called Wilcox-Test)
    wilcox.test(psychic, rock) #data size here too small for EXACT p-value
  
  # Wilcox non-ranked test (non-parametric equivalent of a paired t-test)
    wilcox.test(psychic[1:13], rock, paired = T) # data size here too small for EXACT p-value
  
# linear regression ----
    plot(pokemon$Atk, pokemon$Def)
    regModel <- lm(Def ~ Atk, data = pokemon)   # that's how you plot a linear regression model, has to be supplied with a formula, allows to predict defense 
      regModel                                  # says: we want a linear model lm()
                                                # where Def "depends on" (~) Attack
                                                # and the data is pokemon
                                                # so the Coefficients tell us:
                                                        # intercept: means that when the pokemon's attack = 0, it still has a Defense of 33
                                                        # and for every 1 unit increase in attack --> Defense increases by 0.498
        # this is how you can draw the regression line:
    abline(regModel$coef[1], regModel$coef[2], col = "red", lwd = 3)
                          # coef[1] extracts the intercept
                          # coef[2] extracts the slope
                          # "col" defines color and "lwd" defines line-weight
    
    summary(regModel)     # this is how you get a p-value for the slope
                          # tells you coefficients and the p-values for the coefficients
                          # 1st pr (probability) value tests whether the intercept is = 0
                          # 2nd pr value tests whether the slope = 0 (since it's here very very small it indicates a strong connection b/w ATK/DEF)
                          # multiple R-squared: tells you that 19.8% of the data is accounted for by Atk variable
                          # adjusted R-squared: penalizes you, the more predictors you put into it (the better you try to fit the data in, b/c that would always lead to higher and higher "natural" R-values)
                          # p-Value: tells you whether the model is over all better at predicting Def than the grand "mean"
    hist(regModel$res)  # makes histogram of our regression Model and extracts the "residuals" (res)
    qqnorm(regModel$res)# another way to check normality of your residuals
                        # plots every data point from our data along the sample quantiles and the theoretical quantiles it should follow, IF it was a normal distribution --> should follow a diagonal line
       
    qqline(regModel$res)# make this line visible, use the line function: plots you a line (here red), to show the trend
                        # IF your data was normally distributed, it would follow the created line
      
    # log model might improve Model fit of Defense
    # we will not regress Def depends on Atk, but instead log10Def depends on Atk
    
    regModelLog <- lm(log10(Def) ~ Atk, data = pokemon)  
    regModelLog      
      plot(pokemon$Atk, log10(pokemon$Def))
      abline(regModelLog$coef[1], regModelLog$coef[2], col = "red", lwd = 3) 
          # intercept ~ 1.5 and per Unit Atk, there is a 0.004208 incr. in Def
      summary(regModelLog)
      hist(regModelLog$res, breaks = 40)
        qqnorm(regModelLog$res)
        qqline(regModelLog$res)    
              
      
# ANOVA AND KRUSKAL-WALLIS ----
        
    # we are looking for differences in different types (ghost, grass, ground, and ice)
    # to do that, we are using a subset of the pokemon data
    # makes subset equal to Type Ghost, or Type Grass, or Type Ground, or Type Ice
    pokeSubset        <- pokemon[ pokemon$Type.I == "Ghost" | 
                              pokemon$Type.I == "Grass" |
                              pokemon$Type.I == "Ground"| 
                              pokemon$Type.I == "Ice", ]
                      # looks at the data frame and 
                      # returns every row of these 4 given types
      
        pokeSubset$Type.I <- factor(pokeSubset$Type.I)                      
        plot(pokeSubset[, "Type.I"], pokeSubset[, "Atk"], ylab = "Attack")  
            # to make stuff fit with R (would stop at zeros without this factor fct. --> I guess..)
            # shows the 4 Types + compares Atk, the y-Axis is explicitly called "attack"

    # give me more cum        

        
    # ANOVA oneway Model = linear model to predict Atk based on Type (instead of a continuous predictor, we provide a factor)
              # and the data is pokeSubset
        
      oneWay <- lm(Atk ~ Type.I, data = pokeSubset)   
      oneWay
        # summary(oneWay)
          # can be ignored
          # tests each individual regression, whether ground is diff. from ghost etc.
          # gives multiple/adjusted R-squared
      anova(oneWay)
        
          # ANOVA ~ a special case of regression
          # main point is the p-value (although it also gives you slopes & intercepts & stuff)
            # name of factor: Type I
            # Df, Sum Sq, Mean Sq, F value, Pr <-- P-Value
      
      # --> let's evaluate the model fit, to find out, where the differences to the previous model were!
      
      hist(oneWay$res)
        # looks okay
      qqnorm(oneWay$res)
      qqline(oneWay$res)        
        # pretty good fit
        # assumption: variance between our models is equal --> can be tested with:

  # bartlett test
      bartlett.test(pokeSubset$Atk, pokeSubset$Type.I)
          # we supply the outcome variable --> Atk
          # we supply the grouping variable --> Type

  # One-Way Wallis ANOVA
      oneway.test(Atk ~ Type.I, data = pokeSubset)
        # assumes variance to be unequal: var.equal = FALSE
      oneway.test(Atk ~ Type.I, data = pokeSubset, var.equal = TRUE)
        # assumes var.equal = TRUE --> so we get the same P-value as when we did the lm() fct., followed by the anova fct.


  # Kruskal-Wallis- Rank Sum Test
      # if the data is non-normally distributed, we can use the following test:
        
        kruskal.test(pokeSubset$Atk, pokeSubset$Type.I)
          # gives you the KW-Rank-Sum Test
          # and a P-Value
          
           
# post-hoc tests ----
    TukeyHSD(oneWay)
            # HSD = honestly significant difference
            # compares every group to every other group
            # here: nothing is actually significant..
            # makes sense, since overall p-Value was insignificant as well
        
      # MUCH BETTER: use select pairs of comparisons!
            # first argument:   dependent variable (here Atk)
            # second argument:  predicted variable (here Type I)
   pairwise.t.test(pokeSubset$Atk, pokeSubset$Type.I)     
            # but it's smart to put the p-adjustment method to NONE
            # perorms t-tests between every pairwise comparison and has NOT adjusted the p-values
   pairwise.t.test(pokeSubset$Atk, pokeSubset$Type.I, p.adjust.method = "none")
            # if you're only interested in a specific comparison e.g. grass vs. ghost
            # and ice vs. ground
   p.adjust(c(0.316, 0.069), method = "holm")
            # simply takes the comparisons you're interested in
            # adjusts the p-values for those two comparisons (vs. 6 before)
  

# Two-Way anova ----
   pokeSubset2        <- pokemon[ pokemon$Type.I == "Bug" | 
                                   pokemon$Type.I == "Electric" |
                                   pokemon$Type.I == "Fire"| 
                                   pokemon$Type.I == "Poison", ]
        
   pokeSubset2$Type.I <- factor(pokeSubset2$Type.I)
   boxplot(Atk ~ Captive * Type.I, data = pokeSubset2, col = c("red", 'blue'))
  
   twoway <- lm(Atk ~ Type.I * Captive, data = pokeSubset2)
            # we want to predict, whether our dep.variable Atk is dependent on the Type AND (*) Captive
            # looks at main effect of type, and captive
            # looks at INTERACTION of type AND captive
        twoway
      summary(twoway)    
      anova(twoway)
            
            # now we look if anything changes, if we reverse the order!
    twowayRev <- lm(Atk ~ Captive * Type.I, data = pokeSubset2)
        twowayRev      
      summary(twowayRev)  
      anova(twowayRev)      
            # so the p-values are now different!!!
            # because they are dependent on factors
            # unlike type III ss ( = sums of squares) --> completely independent, whatever order you enter
      
      
# type III ss (sums of squares) --> completely independent ----
      
    twoWay <-     lm(Atk ~ Type.I * Captive, data = pokeSubset2,
                     contrasts = list(Type.I = contr.Sum, Captive = contr.Sum))
    twoWayRev <-  lm(Atk ~ Captive * Type.I, data = pokeSubset2,
                     contrasts = list(Type.I = contr.Sum, Captive = contr.Sum))

    Anova(twoWay, type = 3)  
    Anova(twoWayRev, type = 3)
      # same p-values in both orders!
      # shows complete independence of type III ss
    
    
# post-hoc tests for factorial anova ----
    # stands for post-hoc interaction analysis
    # if you have an interaction, you wanna know if there is a difference for captive and wild for each level of pokemonn type
    
    library(phia)
    testInteractions(twoWay, pairwise = "Captive", fixed = "Type.I")
    
    
    
    
    
    
    
          