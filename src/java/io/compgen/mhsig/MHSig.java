package io.compgen.mhsig;

import io.compgen.cmdline.Help;
import io.compgen.cmdline.License;
import io.compgen.cmdline.MainBuilder;
import io.compgen.common.StringUtils;
import io.compgen.mhsig.cli.SigCall;

import java.io.IOException;

public class MHSig {
	public static String getVersion() {
		try {
			return MainBuilder.readFile("io/compgen/mhsig/VERSION");
		} catch (IOException e1) {
			return "unknown";
		}
	}

	private static String args;
	
	public static String getArgs() {
	    return args;
	}

	public static void main(String[] args) throws Exception {
	    MHSig.args = StringUtils.join(" ", args);
		new MainBuilder()
		.setProgName("mhsig")
		.setHelpHeader("mhsig - Micro-homology signature caller\n---------------------------------------")
		.setDefaultUsage("Usage: mhsig cmd [options]")
		.setHelpFooter("http://compgen.io/mhsig\n"+getVersion())
		.addCommand(Help.class)
		.addCommand(License.class)
		.addCommand(SigCall.class)
		.findAndRun(args);
	}
		
}
