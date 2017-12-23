Refactoring of mtDNA Server code from https://github.com/seppinho/mutation-server such that it takes bam input
from command line.

**Compile**  
javac mtServerCLI.java  
jar cvfe mtServerCLI.jar mtServerCLI *

**Usage**  
java -jar mtServerCLI.jar bam_file

**TO DO**  
Currently assumes rCRS.fasta and rCRS.fasta.fai exist and in the same directory as jar file  
Take command line options  
Other cool stuff  

 



