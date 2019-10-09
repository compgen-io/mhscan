package io.compgen.mhsig.cli;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.IterUtils;
import io.compgen.mhsig.MHSig;
import io.compgen.ngsutils.fasta.IndexedFastaFile;
import io.compgen.ngsutils.vcf.VCFAttributeException;
import io.compgen.ngsutils.vcf.VCFParseException;
import io.compgen.ngsutils.vcf.VCFReader;
import io.compgen.ngsutils.vcf.VCFRecord;

@Command(name="call", desc="Calculate the level of microhomology for indels (from VCF)", category="mhsig")
public class SigCall extends AbstractOutputCommand {
	
	private String genomeFname=null;
	private String vcfFname=null;
	private int minLength = 0;
	private int maxFlankingLength = 1000;
	private boolean passing = false;
	
	private String endKey = null;
	

	@Option(desc="Only use variants that pass all VCF filters", name="passing")
	public void setPassing(boolean val) {
		this.passing = val;
	}

	@Option(desc="Maximum flanking sequence to investigate", name="maxflanking", defaultValue="1000")
	public void setMaxFlankingSeq(int maxFlankingSeq) {
		this.maxFlankingLength = maxFlankingSeq;
	}
	
	@Option(desc="Min indel length", name="minlen", defaultValue="0")
	public void setMinLength(int minLength) {
		this.minLength = minLength;
	}
	
    @Option(desc="Variant file (contains indel coordinates to examine)", name="vcf", required=true)
    public void setVCF(String vcf) {
    	this.vcfFname = vcf;
    }

    @Option(desc="Reference genome FASTA file", name="ref", required=true)
    public void setRef(String ref) {
    	this.genomeFname = ref;
    }

    // TODO: Make this work across BNDs and INVs where you can calculate MH between adjoining sections
    //       Will require alt-chrom and alt-pos to find end points (CHR2, END)
    
    // TODO: Make this work with long range DELs
    //       Will require support for END INFO tag, so add that as an argument

    
    @Option(desc="INFO key for the end of a DEL", name="end-key")
    public void setEndKey(String endKey) {
    	this.endKey = endKey;
    }


	public SigCall() {
	}


	@Exec
	public void exec() throws CommandArgumentException, IOException, VCFParseException {
	    if (genomeFname == null || !(new File(genomeFname)).exists() || !(new File(genomeFname+".fai")).exists()) { 
	        throw new CommandArgumentException("Missing reference FASTA filename (or .fai index)!");
	    }
	    
	    if (vcfFname == null || !(new File(vcfFname)).exists()) { 
	        throw new CommandArgumentException("Missing VCF filename!");
	    }
	    
	    IndexedFastaFile fasta = new IndexedFastaFile(genomeFname);
	    
		VCFReader reader = new VCFReader(vcfFname);

		System.out.println("##mhsig_Command="+MHSig.getArgs());
		System.out.println("##mhsig_Version="+MHSig.getVersion());
		System.out.println("#fasta: "+genomeFname);
		System.out.println("#vcf: "+vcfFname);
		if (passing) {
			System.out.println("#only-passing-variants");
		}
		if (minLength > 0) {
			System.out.println("#min-length: " + minLength);
		}
		
		System.out.println("chrom\tstart\tend\ttype\tref\talt\tlength\tleft_matches\tright_matches\tleft_seq\tright_seq\tindel_seq1\tindel_seq2");
		
		Set<String> invalidAlts = new HashSet<String>();
		invalidAlts.add("<INV>");
		invalidAlts.add("<BND>");
		invalidAlts.add("<DUP>");
		
		for (VCFRecord rec: IterUtils.wrap(reader.iterator())) {
			if (passing && rec.isFiltered()) {
				continue;
			}

			String alt = rec.getAlt().get(0);
			
			if (invalidAlts.contains(alt)) {
				continue;
			}
			
			try {
				if (rec.getInfo().contains("SVTYPE") && rec.getInfo().get("SVTYPE").asString(null).equals("BND")) {
					System.err.println("Skipping BND");
					continue;
				}
			} catch (VCFAttributeException e) {
				e.printStackTrace();
			}
			
			
			String ref = rec.getRef();

			
			int pos = rec.getPos();
			int end = pos;
			
			if (endKey != null && rec.getInfo().contains(endKey)) {
				end = rec.getInfo().get(endKey).asInt();
			}
			
			int delLen = end - pos;
			
			if (rec.getRef().length()>1 || delLen > 1) {
				// this is a deletion
				
				String delLeft="";
				String delRight="";
				
				
				if (delLen == 0) {
					// a normal del, not a delly del...
					
					String del = ref;
					
					while (alt.length()>0 && alt.charAt(0) == del.charAt(0)) {
						alt = alt.substring(1);
						del = del.substring(1);
					}
					
					if (del.length() < minLength) {
						continue;
					}
					
					if (del.length() < maxFlankingLength) {
						delLeft = del;
						delRight = del;
					} else {
						// TODO: fix this... make it only pull subseq (not likley to be used though...)
						delLeft = del;
						delRight = del;
					}
					
					delLen = del.length();
					
				} else {
					
					// for a Delly del, we might need to get the left and right flanks separately.
				
					if (delLen < maxFlankingLength) {
						System.err.println("delLen < maxFlankingLen");
						String del = fasta.fetchSequence(rec.getChrom(), rec.getPos(), rec.getPos()+delLen);
						delLeft = del;
						delRight = del;
					} else {
					
						/*

						====|               |====
						     delL       delR
						
						*/
						 
						
//						System.err.println("delLeft");
						delLeft = fasta.fetchSequence(rec.getChrom(), rec.getPos(), rec.getPos()+maxFlankingLength);
//						System.err.println("delRight");
						delRight = fasta.fetchSequence(rec.getChrom(), rec.getPos()+delLen-maxFlankingLength, rec.getPos()+delLen);
					}
								
				}

				
				/*

                     ====|               |====
				     LLLL                 RRRR
				
				*/

				
				String leftFlank = fasta.fetchSequence(rec.getChrom(), rec.getPos()-delLeft.length(), rec.getPos());
				String rightFlank = fasta.fetchSequence(rec.getChrom(), rec.getPos()+delLen, rec.getPos()+delLen+delRight.length());

				
				
				/*
	     			    LLLL                 RRRR
						====|               |====
						     delL       delR


	     			    LLLL -------------+  RRRR
						====|             V |====
						     delL       delR

				 
	     			    LLLL  +------------- RRRR
						====| V             |====
						     delL       delR

				 */
				
				
				
				// For the left flank, compare the end chars of the leftFlank and delRight.

				int leftMatch = 0;
				for (int i=1; i <= delRight.length(); i++) {
					if (leftFlank.charAt(leftFlank.length()-i) == delRight.charAt(delRight.length()-i)) {
						leftMatch ++;
					} else {
						break;
					}
				}
				
				// For the right flank, compare the start chars of the rightFlank and delLeft.
				int rightMatch = 0;
				for (int i=0; i < delLeft.length(); i++) {
					if (rightFlank.charAt(i) == delLeft.charAt(i)) {
						rightMatch ++;
					} else {
						break;
					}
				}
				
				System.out.println(rec.getChrom()+"\t"+rec.getPos()+"\t"+(rec.getPos()+delLen) + "\tDEL\t" + rec.getRef() + "\t" + rec.getAlt().get(0) + "\t" + delLen + "\t" + leftMatch+"\t" + rightMatch + "\t" + leftFlank + "\t" + rightFlank + "\t" + delLeft + "\t" + delRight);
				
			} else if (alt.length() > 1) {
			
				String ins = alt;
				
				while (ref.length()>0 && ref.charAt(0) == ins.charAt(0)) {
					ref = ref.substring(1);
					ins = ins.substring(1);
				}
				
				if (ins.length() < minLength) {
					continue;
				}
				
				String left = fasta.fetchSequence(rec.getChrom(), rec.getPos()-ins.length(), rec.getPos());
				String right = fasta.fetchSequence(rec.getChrom(), rec.getPos(), rec.getPos()+ins.length());

				int leftMatch = 0;
				
				for (int i=1; i <= ins.length(); i++) {
					if (left.charAt(left.length()-i) == ins.charAt(ins.length()-i)) {
						leftMatch ++;
					} else {
						break;
					}
				}
				
				int rightMatch = 0;
				
				for (int i=0; i < ins.length(); i++) {
					if (right.charAt(i) == ins.charAt(i)) {
						rightMatch ++;
					} else {
						break;
					}
				}
				
				System.out.println(rec.getChrom()+"\t"+rec.getPos()+"\t"+rec.getPos() + "\tINS\t" + rec.getRef() + "\t" + rec.getAlt().get(0) + "\t" + ins.length() + "\t" + leftMatch+"\t" + rightMatch + "\t" + left + "\t" + right + "\t" + ins);
				
			}
		}		

		reader.close();
	    
	}
}
