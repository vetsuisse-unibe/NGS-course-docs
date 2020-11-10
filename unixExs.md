## Exercises 
In this Practical we will login into the Bioinformatics server (binfservms01.unibe.ch) and do all the exercises on the server.  

In order to login into the Bioinformatics cluser you need to be connected to unibe network via VPN or eduroam. 
### VPN 
If you are working from home then please login into the unibe VPN server using cisco Anyconnect or FortiConnect. 
(https://tinyurl.com/yxj8ovpd)

### eduroam
If you are connected to the eduroam network then you don't need to connect to VPN server.  

### SSH into bionformatics server 
One you are connected to the VPN server or eduroam network, ssh into server in the follow manner according to your OS.  
Check the Login info word document for login and password details.  

### Mac OS X users 
Locate your Terminal.app and type the following *ssh command*
```
ssh login@binfservms01.unibe.ch
```
When it prompts for password please enter the corresponding password also from the LoginInfo document.
### Windows 10 Users 
Locate Ubuntu and type the following *ssh command*
```
ssh login@binfservms01.unibe.ch
```
When it prompts for password please enter the corresponding password also from the LoginInfo document.

### Older Windows Users 
Locate Mobaxterm and type the following *ssh command* 
```
ssh login@binfservms01.unibe.ch
```
When it prompts for password please enter the corresponding password also from the LoginInfo document.

#### Please ask for assistance any time, you do not understand the exercises. 

Commands in grey blocks are to be typed in the Unix command shell or command line prompt.  

Words in italics need to be replaced by the proper parameters (for example, your file name).

### Try some basic unix commands
*Display user name*
```
whoami 
```
*Show the current working directory*
```
pwd 
```
It should show something like /home/student41. Which is 'home' directory for user student01

*List the files in the directory*	
```
ls 
```
*Create a empty file*	
```
touch <filename> 
```
*Print a string to the screen*
```
echo "Hello world" 
```
*Print the current date*
```
date
```
*See a history of all the last commands you tried* 
```
history 
```
*Get local help page of a command* 
```
man ls
```
*Run the following commands one after another*
```
touch <exampleFile>
ls
```
*Print the contents of the file* 
```
cat <exampleFile>
```
It should show no lines as it is an empty file. 

*Rename a file* 
```
mv <exampleFile> <exampleFile2>
```
*Delete a file* 
```
rm <exampleFile2>
```
*Create a new folder/directory* 
```
mkdir <exampleDirectory> 
```
*Create a file under the new folder/directory*
```
touch exampleDirectory/exampleFile
```
*Delete the folder/directory* 
```
rmdir exampleDirectory
```
worked ? No ! 

*Delete the file first* 
```
rm exampleDirectory/exampleFile
```
*Delete the directory now* 
```
rmdir exampleDirectory
```
# creating and moving around Directories 
*Create a directory called 'Documents' and change current directory to Documents* 
```
mkdir Documents 
cd Documents 
pwd
```
*Going up one directory* 
```
cd..
```
*and then type to see what has happened*
```
pwd
```
*Go up by two directories* 
```
cd ../.. 
```
*Go to home directory* 
```
cd
```
Always type “pwd” to locate yourself

# Command arguments

Most programs in UNIX accept arguments that modify the program’s behavior. For example 
List the files in longer format 
```
ls -l 
touch exampleFile1
touch exampleFile2

ls
ls -l 
```
Different example parameters used with ls 
```
ls -a List all files, including hidden ones.
ls -h List all files, with human-readable sizes (Mb, Gb).
ls -l List all files, long format.
ls -S List all files, order by size.
ls -t List all files, order by modification time.
ls -1 List all files, one file per line.
```

Parameters for remove command 
```
rm exampleFile1
rm -i exampleFile2
```
Aliases are short forms used for commands. 
```
alias rm=”rm -I”
touch <exampleFile>
rm <exampleFile>
```
So better to alias rm as rm -i to be on the safer side ?.


# Redirection.

All the above commands sent the output if any to the screen. Instead of outputting on the screen redirection helps you to put into a file 
```
echo "My first line" >testFile.txt
cat testFile.txt 
echo "My second line" >>testFile.txt
cat testFile.txt
ls / >> ListRootDir.txt
cat ListRootDir.txt
```

# Wildcards

wildcard is a symbol that is used to represent one or more characters. Example wildcards are as follows 
* Zero or more characters
? Any single character. 

In the following examples you can test some of these wildcard characters 

Create a new folder 
mkdir wildCardTesting

Change directory 
```
cd wildCardTesting
touch test1.txt 
touch test2.txt 
```
Create several files with a single command 
```
touch test3.txt test4.txt test1.csv test2.csv test3.csv test4.csv
```
Count the number of files using the pipes ( | symbol) 
wc –l counts the number of lines in the input  
```
ls -l | wc -l 
```

Try using the following wild cards 
```
ls * 
ls test*
ls *.txt 
ls *.csv
ls test[1-2].txt 
ls test[!3].*
```

# subset a file with grep and awk
```
mkdir GenomeStats 
cd GenomeStats
```
Download the text file showing the submitted Genomes for different Eukaryote species at NCBI https://www.ncbi.nlm.nih.gov/genome/browse/
```
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"
```
Use less to have a quick view of the file 
```
less eukaryotes.txt
```
It is a tab delimited text file with several columns. 
The first line shows the different column headers
```
head -n 1 eukaryotes.txt
```
We want to see how many cow assemblies have been submitted 
```
grep "Bos taurus" eukaryotes.txt
```
More easier ^ stands for beginning of  a line 
```
grep "^Bos taurus" eukaryotes.txt | wc -l 
grep -c "^Bos taurus" eukaryotes.txt
```
Want to keep the header line ? 
```
grep -E '^#|^Bos taurus' eukaryotes.txt
```
Lets do some statistics on available Genomes 
How many Animal and plant genomes are available 
cut command in unix can be used to select columns from tab de-limited files 
Only the column Group can selected using cut 
```
cut -f5 eukaryotes.txt | less 
```
Now the pipes can be used to see the number of animal plant genomes available at NCBI 
```
cut -f 5 eukaryotes.txt |sort | uniq -c
```
So how many animal and plant genomes 
Now use cut and pipe symbol to find the number of mammalian genomes available at NCBI. (Hint: check column 6) 
cut can be used to select more columns 
```
cut -f 1,6,8 eukaryotes.txt |less
```
which Mammalian genome has the highest GC content 
```
cut -f 1,6,8 eukaryotes.txt | grep "Mammals" | sort -t$'\t' -nrk3 |less
cut -f 1,6,8 eukaryotes.txt | grep "Mammals" | sort -t$'\t' -nrk3 | head –n 1
```
which Mammalian genome has the least GC content 

# Question of the day
Is the statement "a genome-wide GC content of ≈30% is one of the lowest observed in any animal genome"  True ? 

If the find the answer for this you just proved or disproved an accepted hypothesis!

# Welcome to the exciting world of Data Analysis.