# SOME BUILT-IN FUNCTIONS ----
normalDist <- rnorm(1000, 0, 1)
mean(normalDist)
hist(normalDist)
?hist
hist(normalDist, breaks = 50)
 
# WRITING HELLO FUNCTION ----
sayHello <- function(){
  "Hello!"
}

sayHello()

sayHello <- function(name){
  paste("Hello", name)
}

greet <- sayHello("Hefin")
greet

# FUNCTION WITH OPTIONAL ARGUMENTS ----
expo <- function(x, exp = 2){
  hist(x ^ exp, ...)
}

expo()
expo(normalDist)
expo(normalDist, 1, breaks = 2)

# FUNCTION WITH UNNAMED ARGUMENTS ----
expo <- function(x, exp = 2, ...){
  hist(x ^ exp, ...)
}

expo()
expo(normalDist)
expo(normalDist, 1, breaks = 2)

# FUNCTION WITH LOGICAL ARGUMENTS ----
expo <- function(x, exp = 2, hist = FALSE, ...){
  if(hist == TRUE){
    hist(x ^ exp, ...)
    x ^ exp
  } else{
    x ^ exp
  }
}

expo(normalDist)
expo(normalDist, hist = TRUE)

# POPULATION DOUBLING COUNTER ----
## 10,000 cells -> 20,000 cells = 1 population double
## 10,000 cells -> 40,000 cells = 2 population doubles
popDouble <- function(start, end){
  doubling <- -1
  while(end > start) {
    doubling <- doubling + 1
    start <- 2 * start
  }
  doubling + (end / start)
}

popDouble(10000, 20000)
popDouble(10000, 30000)
popDouble(10000, 40000)
popDouble(10000, 1000000)
