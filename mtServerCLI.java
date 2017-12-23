//package genepi.mut.pileup;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;

//import genepi.base.Tool;
//import genepi.io.FileUtil;
import genepi.io.text.LineWriter;

import genepi.mut.objects.BasePosition;
import genepi.mut.objects.VariantLine;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndex;
//import htsjdk.samtools.reference.FastaSequenceIndexCreator;



public class mtServerCLI{

	String version = "mtdna";
	//String input;

	File bam_file;
	String refPath;

	String file_sample_id;
	String outputRaw;
	String outputVar;

	String indel;
	String baq;

	double level;
	int baseQ;
	int mapQ;
	int alignQ;

	public mtServerCLI(String filename){
		refPath = "rCRS.fasta";
		indel = "true";
		baq = "false";

		level = 0.01;
		baseQ = 20;
		mapQ = 20;
		alignQ = 30;

		bam_file = new File(filename);
	    File ref_file = new File(refPath);
	    File ref_file_index = new File(refPath+".fai");

	    if(!bam_file.exists()){
	      System.out.println("Bam file: "+filename+ " does not exits!");      
	      System.exit(1);
	    }

	    if(!ref_file.exists() || !ref_file_index.exists()){
	      System.out.println("Ref file "+refPath+" or "+refPath+".fai does not exits!");
	      System.exit(1);                  
	    }

    	SAMFileReader inputSam = new SAMFileReader(bam_file);
    	String sample_id = inputSam.getFileHeader().getReadGroups().get(0).getSample();

    	file_sample_id = sample_id.replaceAll(" ", "_");    	
    	outputRaw = file_sample_id+"_raw.txt";
    	outputVar = file_sample_id+"_var.txt";
    	
	}

	public int run() {
		
		LineWriter writerRaw = null;
		LineWriter writerVar = null;

		try {			
			File outRaw = new File(outputRaw);
			File outVar = new File(outputVar);

			writerRaw = new LineWriter(outRaw.getAbsolutePath());
			writerRaw.write(BamAnalyser.headerRaw);

			writerVar = new LineWriter(outVar.getAbsolutePath());
			writerVar.write(BamAnalyser.headerVariants);

		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		long start = System.currentTimeMillis();

		/*
		BamAnalyser analyser = new BamAnalyser(bam_file.getName(), refPath, baseQ, mapQ, alignQ, Boolean.valueOf(baq),
				version);
		*/
		BamAnalyser analyser = new BamAnalyser(file_sample_id, refPath, baseQ, mapQ, alignQ, Boolean.valueOf(baq),
				version);

		System.out.println("Processing: " + bam_file.getName());

		try {
			analyseReads(bam_file, analyser, Boolean.valueOf(indel));
			determineVariants(analyser, writerRaw, writerVar, level);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return 0;
		}

		try {
			writerVar.close();
			writerRaw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println("Raw file written to " + new File(outputRaw).getAbsolutePath());
		System.out.println("Variants file written to " + new File(outputVar).getAbsolutePath());
		System.out.println("Time: " + (System.currentTimeMillis() - start) / 1000 + " sec");
		return 0;
	}
	

	
	// mapper
	private void analyseReads(File file, BamAnalyser analyser, boolean indelCalling) throws Exception, IOException {

		// TODO double check if primary and secondary alignment is used for
		// CNV-Server
		final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT)
				.open(file);

		SAMRecordIterator fileIterator = reader.iterator();

		while (fileIterator.hasNext()) {
			SAMRecord record = fileIterator.next();
			analyser.analyseRead(record, indelCalling);
		}
		reader.close();
	}
	

	// reducer
	private void determineVariants(BamAnalyser analyser, LineWriter writerRaw, LineWriter writerVar, double level)
			throws IOException {

		HashMap<String, BasePosition> counts = analyser.getCounts();

		String reference = analyser.getReferenceString();

		for (String key : counts.keySet()) {
			String idKey = key.split(":")[0];
			String positionKey = key.split(":")[1];
			int pos;
			boolean insertion = false;

			if (positionKey.contains(".")) {
				pos = Integer.valueOf(positionKey.split("\\.")[0]);
				insertion = true;
			}
			else{
				pos = Integer.valueOf(positionKey);
			}

			if (pos > 0 && pos <= reference.length()) {
				char ref = 'N';
				BasePosition basePos = counts.get(key);
				basePos.setId(idKey);
				basePos.setPos(pos);

				VariantLine line = new VariantLine();

				if (!insertion) {
					ref = reference.charAt(pos - 1);

				} 
				else{
					line.setInsPosition(positionKey);
				}

				line.setRef(ref);
				line.analysePosition(basePos, level);
				line.callVariants(level);

				if (line.isFinalVariant()) {
					writerVar.write(line.writeVariant());
				}

				// raw data
				String raw = line.toRawString();
				writerRaw.write(raw);

			}

		}
	}

	public static void main(String[] args) {

		if(args.length != 1){
	      System.out.println("Usage: java -jar mtServerCLI.jar bam_file");
	      System.exit(1);
	    }

	    mtServerCLI pileup = new mtServerCLI(args[0]);
	    pileup.run();



	}

}
