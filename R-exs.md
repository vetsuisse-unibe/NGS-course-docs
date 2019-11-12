## R Exercises 

Let's work with R using RStudio. 

So start RStudio in your computer. Once its open you will see the 4 panels in RStudio window: 
* the text editor for your scripts and documents (top-left, in the default layout), 
* the R Console (bottom-left), 
* your Environment/History (top-right),  
* your Files/Plots/Packages/Help/Viewer (bottom-right). 

The following codes should be typed in the console. 

#### Creating objects in R 
#You can get output from R simply by typing math in the console 
```{r, eval = FALSE}
16 + 16 
14/7
(2^3 - 6) / 2
```
#To do more interesting things it is better to assin values to objects as follows. 

```
x <-55
```
#<- is the assignment, the value is on the right and is assigned or stored in the named object on the left. 
#Now the object x is in memory and any kind of arithmetic can be done 
```
2 * x 
2 + x 
2 - x 
```
#The arithmetic results can be assigned to a new object or to the same object 
```
y <-2*x 
y 
x<-2*x 
x 
```
#Create a vector using the function c()
```
animals <-c("cat","dog","cow")
genomeSizesMB <-c(2641.34, 2410.98,2670.14)
```
Note Character vector has "" quotes and the numerical vector has none. 

Data Frame is another object commonly used in R for tabular data. It can be seen as a collection of vectors. The vectors can be only character vectors or only number vectors or a mix of character and numerical vectors. 
```
animals <-c("cat","dog","cow")
genomeSizesMB <-c(2641.34, 2410.98,2670.14)
genes <-c(34949, 36809, 32432)
centre<-c("International Cat Genome Sequencing Consortium", "Dog Genome Sequencing Consortium", "Center for Bioinformatics and Computational Biology, University of Maryland") 
df <-data.frame(animals, genomeSizesMB,genes,centre)
df
str(df) 
```
The output of str shows that animals and centre are classified as **factor vectors** 

you can also use View to check the data in the GUI fashion

```
View(df)
```
matrix is another object in R. Matrix can be called as a numerical data frame 
```
n <- c(1:9)
mat <-matrix(n,nrow=3,ncol=3,byrow=TRUE)
mat
```
#### Subsetting with objects 
if we want extract one or two values from vectors it can be done with vector name followed by square brackets 

```
animals <-c("cat","dog","cow")
animals[1]
animals[2]
animals[3]
````
For extracting more values it can be done as follows 
```
animals[c(1,3)]
```
To retrieve data from matrix or dataframe we would enter its row and column coordinates in the single square bracket "[]" operator. Lets extract from the df dataframe that you created before 

extracts the genomeSize of cat in MB
```
df[1,2] 
```
extracts the first column which are animals
```
df[,1] 
```
extracts the first row
```
df[1,] 
```

#### Functions and arguments 
Many statistical methods like mean, standard deviation etc are already made available in R which can be used readily 

```
a <-49 
b<-sqrt(49)
mean(genes)
```
#### R Packages 
In the following lines you will install a R package called dplyr. Dplyr is R package used to transform and summarize tabular data. 

```
install.packages("dplyr")
# load the package
library(dplyr)
data<-read.delim2("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS//eukaryotes.txt",sep="\t",header=T)
str(data)
```

The following code converts some of the columns read as factors into numeric 
```
data$GC <-as.numeric(as.character(data$GC))
data$Size..Mb.<-as.numeric(as.character(data$Size..Mb.))
data$TaxID <-as.numeric(as.character(data$TaxID))
data$Genes<-as.numeric(as.character(data$Genes))
data$Proteins <-as.numeric(as.character(data$Proteins))
```
filter by rows only Land Plants data 
```
data.plants <- filter(data,SubGroup=="Land Plants")
summary(data.plants) 
```

Visualising the number of genes in Land Plants 
```
hist(data.plants$Genes, xlab="Land Plants Proteins") 
```
#### GC Graphs 
Use the eukaryotes.txt that you downloaded previously to make a comparison graph of the %GC in different SubGroups
```
library(lattice)
densityplot(~data$GC, groups=data$SubGroup,data=data,auto.key=TRUE)
#you can add the legend below to make the plot look better 
densityplot(~data$GC, groups=data$SubGroup,data=data,auto.key=list(space="bottom", columns=2, title="SubGroup",cex.title=0.6, cex=.5))
densityplot(~data$GC|data$SubGroup,data=data)
```
Question : Which subgroup has more varied %GC. 






