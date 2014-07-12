package dnatoproteins;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import javax.swing.JOptionPane;

import org.apache.commons.io.FileUtils;
import org.biojava3.core.sequence.RNASequence;

public class DNAtoProtein {
	  // List of characters being used
	  private static final String CHAR_LIST = 
		        "ATCG";
	        // Maximum length of DNA generated in the RandomDNA.txt file
		    private static final int RANDOM_STRING_LENGTH = 103;
		    
		    // Randomly generate either a '+' or a '-' to determine transcription direction
		    public String randomSymbol(){
		    	ArrayList<String> symbols = new ArrayList<String>();
		    	symbols.add("+");
		    	symbols.add("-");
		    	Random random = new Random();
		    	return symbols.get(random.nextInt(2)).toString();
		    }
		     
		    /**
		     * This method generates random DNA strand
		     * @return
		     */
		    public String generateRandomString(){
		         
		        StringBuffer randStr = new StringBuffer();
		        for(int i=0; i<RANDOM_STRING_LENGTH; i++){
		            int number = getRandomNumber();
		            char ch = CHAR_LIST.charAt(number);
		            randStr.append(ch);
		        }
		        return randStr.toString();
		    }
		     
		    /**
		     * This method generates random numbers
		     * @return int
		     */
		    private int getRandomNumber() {
		        int randomInt = 0;
		        Random randomGenerator = new Random();
		        randomInt = randomGenerator.nextInt(CHAR_LIST.length());
		        if (randomInt - 1 == -1) {
		            return randomInt;
		        } else {
		            return randomInt - 1;
		        }
		    }
		    
		    // This method generates the corresponding RNA strand depending on the String value passed in the constructor
		    public String generatemRNA(String dna){
	    		StringBuilder sb = new StringBuilder();
		    	if(dna.endsWith("+")){
		    		char[] chars = dna.toCharArray();
		    		for(int x = 0; x < chars.length; x++){
		    			if(chars[x] == 'A'){
		    				chars[x] = 'U';
		    				sb.append(chars[x]);
		    			}
		    			else if(chars[x] == 'T'){
		    				chars[x] = 'A';
		    				sb.append(chars[x]);
		    			}
		    			else if(chars[x] == 'C'){
		    				chars[x] = 'G';
		    				sb.append(chars[x]);
		    			}
		    			else if(chars[x] == 'G'){
		    				chars[x] = 'C';
		    				sb.append(chars[x]);
		    			}
		    		}
		    	}
		    	if(dna.endsWith("-")){
		    		char[] chars = dna.toCharArray();
		    		for(int x = 0; x < chars.length; x++){
		    			if(chars[x] == 'A'){
		    				chars[x] = 'U';
		    				sb.append(chars[x]);
		    			}
		    			else if(chars[x] == 'T'){
		    				chars[x] = 'A';
		    				sb.append(chars[x]);
		    			}
		    			else if(chars[x] == 'C'){
		    				chars[x] = 'G';
		    				sb.append(chars[x]);
		    			}
		    			else if(chars[x] == 'G'){
		    				chars[x] = 'C';
		    				sb.append(chars[x]);
		    			}
		    		}
		    		String prerna = sb.reverse().toString();
		    		return prerna;
		    	}
		    	return sb.toString();
		    }
		    // This method splits the string every n times. Used for getting codons.
		    public String[] splitStringEvery(String s, int interval) {
		        int arrayLength = (int) Math.ceil(((s.length() / (double)interval)));
		        String[] result = new String[arrayLength];

		        int j = 0;
		        int lastIndex = result.length - 1;
		        for (int i = 0; i < lastIndex; i++) {
		            result[i] = s.substring(j, j + interval);
		            j += interval;
		        } //Add the last bit
		        result[lastIndex] = s.substring(j);

		        return result;
		    }
		    
		    private static final Map<String ,String> TRANSLATION = new HashMap<String ,String>(); 
		    
		    // Codons to translate 
		    private static final String[] CODONS = { "uuu", "uuc", "uua", "uug", 
		    "ucu", "ucc", "uca", "ucg", 
		    "uau", "uac", "uaa", "uag", 
		    "ugu", "ugc", "uga", "ugg", 
		    "cuu", "cuc", "cua", "cug", 
		    "ccu", "ccc", "cca", "ccg", 
		    "cau", "cac", "caa", "cag", 
		    "cgu", "cgc", "cga", "cgg", 
		    "auu", "auc", "aua", "aug", 
		    "acu", "acc", "aca", "acg", 
		    "aau", "aac", "aaa", "aag", 
		    "agu", "agc", "aga", "agg", 
		    "guu", "guc", "gua", "gug", 
		    "gcu", "gcc", "gca", "gcg", 
		    "gau", "gac", "gaa", "gag", 
		    "ggu", "ggc", "gga", "ggg" 
		    }; 
		    
		    // Amino acid in map 
		    private static final String[] AMINO_ACIDS = { "F", "F", "L", "L", 
		    "S", "S", "S", "S", 
		    "Y", "Y", "--STOP--", "--STOP--", 
		    "C", "C", "--STOP--", "W", 
		    "L", "L", "L", "L", 
		    "P", "P", "P", "P", 
		    "H", "H", "Q", "Q", 
		    "R", "R", "R", "R", 
		    "I", "I", "I", "M", 
		    "T", "T", "T", "T", 
		    "N", "N", "K", "K", 
		    "S", "S", "R", "R", 
		    "V", "V", "V", "V", 
		    "A", "A", "A", "A", 
		    "D", "D", "E", "E", 
		    "G", "G", "G", "G" 
		    };
		    
		    public static void init() { 
		    	 for (int i=0; i<CODONS.length; i++) 
		    	 TRANSLATION.put(CODONS[i], AMINO_ACIDS[i]); 
		    }
		    
		    public String generateProteinSequence(String line){
		    	String seq = new RNASequence(line).getProteinSequence().toString();
		        return seq;
		    }


	public static void main(String[] args) throws IOException {		
		
		JOptionPane.showMessageDialog(null, "Click ok to begin");
		
		DNAtoProtein main = new DNAtoProtein();
        String path = "C:"+File.separator+"DNAtoProtein"+File.separator+"RandomDNA.txt";
        final File inputfile = new File(path);
        // creating the file
        inputfile.createNewFile();
        
        // writing the randomly generated DNA strands to the RandomDNA.txt file
        FileWriter writer = new FileWriter(inputfile.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(writer);
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(",	"+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(", "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(" "+main.randomSymbol());
        bw.newLine();
        bw.write(main.generateRandomString());
        bw.write(" "+main.randomSymbol());
        bw.close();
        
        
        // Successful Operation
        JOptionPane.showMessageDialog(null, "Random DNA strands generated successfully. Click Ok to continue.");
        JOptionPane.showMessageDialog(null, "Click OK to begin generating the corresponding mRNA strands");
       // Creating the path for the mRNA.txt file
       String outpath = "C:"+File.separator+"DNAtoProtein"+File.separator+"mRNA.txt";

       
       File output = new File(outpath);
       // creating the mRNA.txt file
       output.createNewFile();
       
       // Writing the mRNA strands into the mRNA.txt file
       FileWriter fw = new FileWriter(output.getAbsoluteFile());
       BufferedWriter rna = new BufferedWriter(fw); 
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(0)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(1)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(2)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(3)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(4)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(5)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(6)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(7)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(8)));
       rna.newLine();
       rna.write(main.generatemRNA(FileUtils.readLines(inputfile).get(9)));
       rna.close();
       
       // Successful Operation
       JOptionPane.showMessageDialog(null, "DNA transcribed successfully!");
       JOptionPane.showMessageDialog(null, "Click OK to begin mRNA translation to proteins");
       // Creating the path for the ProteinSequence.txt file
       String proteinpath = "C:"+File.separator+"DNAtoProtein"+File.separator+"ProteinSequence.txt";
       
       File proteins = new File(proteinpath);
       // Creating the ProteinSequence.txt file
       proteins.createNewFile();
       
       // And, finally, writing the protein sequences on the ProteinSequence.txt file
       FileWriter lefw = new FileWriter(proteins.getAbsoluteFile());
       BufferedWriter proteinwriter = new BufferedWriter(lefw); 
       proteinwriter.write("***NOTE*** A '*' represents a STOP codon!");
       proteinwriter.newLine();
       proteinwriter.write("It's open source! Created a github repository. Link: https://github.com/MoeMixMC/DNA-to-Protein-Sequencer");
       proteinwriter.newLine();
       proteinwriter.write("Protein Sequences:");
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(0)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(1)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(2)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(3)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(4)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(5)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(6)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(7)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(8)));
       proteinwriter.newLine();
       proteinwriter.write(main.generateProteinSequence(FileUtils.readLines(output).get(9)));
       proteinwriter.close();
       
       // Last successful operation notification :)
       JOptionPane.showMessageDialog(null, "Protein Sequences generated succesfully!");
	}

}
