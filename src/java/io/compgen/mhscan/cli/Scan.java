package io.compgen.mhscan.cli;

import java.io.File;
import java.io.IOException;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;
import io.compgen.common.IterUtils;
import io.compgen.mhscan.MHScan;
import io.compgen.ngsutils.fasta.IndexedFastaFile;
import io.compgen.ngsutils.vcf.VCFParseException;
import io.compgen.ngsutils.vcf.VCFReader;
import io.compgen.ngsutils.vcf.VCFRecord;

@Command(name="scan", desc="Calculate the level of microhomology for indels (from VCF)", category="mhscan")
public class Scan extends AbstractOutputCommand {
	
	private String genomeFname=null;
	private String vcfFname=null;
	private int minLength = 0;
	private boolean passing = false;

	@Option(desc="Only use variants that pass all VCF filters", name="passing")
	public void setPassing(boolean val) {
		this.passing = val;
	}

	@Option(desc="Min indel length", name="minlen", defaultValue="0")
	public void setMinLength(int minLength) {
		this.minLength = minLength;
	}
	
    @Option(desc="Variant file (contains indel coordinates to examine)", name="vcf")
    public void setVCF(String vcf) {
    	this.vcfFname = vcf;
    }

    @Option(desc="Reference genome FASTA file", name="ref")
    public void setRef(String ref) {
    	this.genomeFname = ref;
    }
    
	public Scan() {
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

		System.out.println("##mhscan_Command="+MHScan.getArgs());
		System.out.println("##mhscan_Version="+MHScan.getVersion());
		System.out.println("#fasta: "+genomeFname);
		System.out.println("#vcf: "+vcfFname);
		if (passing) {
			System.out.println("#only-passing-variants");
		}
		if (minLength > 0) {
			System.out.println("#min-length: " + minLength);
		}
		
		System.out.println("chrom\tstart\tend\ttype\tref\talt\tlength\tleft_matches\tright_matches\tleft_seq\tindel_seq\tright_seq");
		
		for (VCFRecord rec: IterUtils.wrap(reader.iterator())) {
			if (passing && rec.isFiltered()) {
				continue;
			}

			String alt = rec.getAlt().get(0);
			String ref = rec.getRef();

			if (rec.getRef().length()>1) {
				// this is a deletion
				
				String del = ref;
				
				while (alt.length()>0 && alt.charAt(0) == del.charAt(0)) {
					alt = alt.substring(1);
					del = del.substring(1);
				}
				
				if (del.length() < minLength) {
					continue;
				}
				
				String left = fasta.fetchSequence(rec.getChrom(), rec.getPos()-del.length(), rec.getPos());
				String right = fasta.fetchSequence(rec.getChrom(), rec.getPos()+del.length(), rec.getPos()+del.length()+del.length());

				int leftMatch = 0;
				
				for (int i=1; i <= del.length(); i++) {
					if (left.charAt(left.length()-i) == del.charAt(del.length()-i)) {
						leftMatch ++;
					} else {
						break;
					}
				}
				
				int rightMatch = 0;
				
				for (int i=0; i < del.length(); i++) {
					if (right.charAt(i) == del.charAt(i)) {
						rightMatch ++;
					} else {
						break;
					}
				}
				
				
				System.out.println(rec.getChrom()+"\t"+rec.getPos()+"\t"+(rec.getPos()+del.length()) + "\tDEL\t" + rec.getRef() + "\t" + rec.getAlt().get(0) + "\t" + del.length() + "\t" + leftMatch+"\t" + rightMatch + "\t" + left + "\t" + del + "\t" + right);
				
				
				
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

				
				System.out.println(rec.getChrom()+"\t"+rec.getPos()+"\t"+rec.getPos() + "\tINS\t" + rec.getRef() + "\t" + rec.getAlt().get(0) + "\t" + ins.length() + "\t" + leftMatch+"\t" + rightMatch + "\t" + left + "\t" + ins + "\t" + right);

				
			}
		}		

		reader.close();
	    
	}
}
