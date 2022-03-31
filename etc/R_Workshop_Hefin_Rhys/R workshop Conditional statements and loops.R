# IFELSE ----
ifelse(4 > 5, "Yes!", "No!")
ifelse(5 >= 4, "Yes!", "No!")

# IF STATEMENTS ----
a <- 4
b <- 5
if(a < b){
  b / a
}

# ELSE STATEMENTS ----
if(a < b){
  paste("b / a = ", b / a)
} else{
    paste("a / b = ", a / b)
}

# ELSEIF STATEMENTS ----
if(a < b){
  paste("b / a = ", b / a)
 } else if(a == b){
    "a and b are equal"
} else{
    paste("a / b = ", a / b)
}

# WHILE STATEMENTS ----
while(1 == 1){
  print("This is true")
}

x <- 1
y <- 10

while(x < y){
  print("This is true")
  x <- x + 1
}

# FOR LOOPS ---
result <- c()

for(i in 1:5){
  result[i] <- i ^ 2
}

result

# COMBINING IF AND FOR ----
for(i in 1:5){
  if(i %% 2 == 0){
    print("Even")
  } else{
    print("Odd")
  }
}

# LAPPLY ----
myList <- list()
for(i in 1:1000){
  norm <- list(rnorm(10, 0, 1))
  myList <- c(myList, norm)
}

head(myList, 3)

listMeans <- list()

for(sample in 1:length(myList)){
  listMeans[sample] <- mean(myList[[sample]])
}

head(listMeans)

listMeans2 <- lapply(myList, mean)
head(listMeans2)

# SAPPLY ----
listMeans3 <- sapply(myList, mean)
listMeans3

# APPLY ----
data(iris)
apply(iris[, 1:4], 2, mean)
apply(iris[, 1:4], 1, mean)
