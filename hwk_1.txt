Lecture 1
1. Done 
2. Done 
3. After installing a software, i would get a message "Biostar Handbook Installation Complete" and that is how i knew the installation was successful 
4. Yes
5. Version: 1.20 (using htslib 1.20)
6. https://github.com/niralds/BMMB-852

Lecture 2
1. a unix command i would find helpful is the open command. So far we looked at how to view list of files in a directory, create directories, rename them and create a txt file. But we didn't learn how to open the file. With terminal, i can easily open the files i need to use without having to clicking buttons. 
2.using the flag command, where i type in "open -p ~/work/test/opening_file.txt will give me an error message stating "invalid option". But if i use "open -e ~/work/test/opening_file.txt" the file will open. I found this by searching the manual "man open"
3. -s flag will sort the files from largest to smallest, and -l gives the long format of a file and the size of the file would be written. When you read through it, it is clear which number refers to the file size 
4. -i flag will request for confirmation before attempting to remove each file
5. nested file of creating edu directory in work and then creating a txt file in edu directory
nshah@Nirals-MacBook-Air ~
$ mkdir -p ~/work/edu
nshah@Nirals-MacBook-Air ~
$ cd work/edu
nshah@Nirals-MacBook-Air ~/work/edu
$ touch bio_1.txt
nshah@Nirals-MacBook-Air ~/work/edu
$ ls
bio_1.txt
nshah@Nirals-MacBook-Air ~/work/edu
$

6. absolute and relative path to a file

nshah@Nirals-MacBook-Air ~/work/edu
$ cd ../test/niral/
nshah@Nirals-MacBook-Air ~/work/test/niral
$ cd
nshah@Nirals-MacBook-Air ~
$ cd work/test/niral/
nshah@Nirals-MacBook-Air ~/work/test/niral
$
demonstrating directory shortcuts
 nshah@Nirals-MacBook-Air ~/work/test/niral
$ cd ./
nshah@Nirals-MacBook-Air ~/work/test/niral
$ cd ../
nshah@Nirals-MacBook-Air ~/work/test
$ cd ~/
nshah@Nirals-MacBook-Air ~
$