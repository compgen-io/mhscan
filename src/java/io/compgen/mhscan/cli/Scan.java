package io.compgen.mhscan.cli;

import io.compgen.cmdline.annotation.Command;
import io.compgen.cmdline.annotation.Exec;
import io.compgen.cmdline.annotation.Option;
import io.compgen.cmdline.annotation.UnnamedArg;
import io.compgen.cmdline.exceptions.CommandArgumentException;
import io.compgen.cmdline.impl.AbstractOutputCommand;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

@Command(name="scan", desc="Calculate the level of microhomology for indels (from VCF)", category="mhscan")
public class Scan extends AbstractOutputCommand {
	
	private String genomeFname=null;
	private String vcfFname=null;
	
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
	public void exec() throws CommandArgumentException, IOException {
    if (genomeFname == null) { 
        throw new CommandArgumentException("Missing --ref argument!");
    }
    if (vcfFname == null) { 
        throw new CommandArgumentException("Missing --vcf argument!");
    }
 }
}
