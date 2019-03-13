package io.compgen.mhscan;

import io.compgen.mhscan.cli.Scan;
import io.compgen.cmdline.Help;
import io.compgen.cmdline.License;
import io.compgen.cmdline.MainBuilder;
import io.compgen.common.StringUtils;

import java.io.IOException;

public class MHScan {
	public static String getVersion() {
		try {
			return MainBuilder.readFile("io/compgen/mhscan/VERSION");
		} catch (IOException e1) {
			return "unknown";
		}
	}

	private static String args;
	
	public static String getArgs() {
	    return args;
	}

	public static void main(String[] args) throws Exception {
	    MHScan.args = StringUtils.join(" ", args);
		new MainBuilder()
		.setProgName("mhscan")
		.setHelpHeader("mhscan - Micro-homology scanner\n---------------------------------------")
		.setDefaultUsage("Usage: mhscan cmd [options]")
		.setHelpFooter("http://compgen.io/mhscan\n"+getVersion())
		.addCommand(Help.class)
		.addCommand(License.class)
		.addCommand(Scan.class)
		.findAndRun(args);
	}
		
}
