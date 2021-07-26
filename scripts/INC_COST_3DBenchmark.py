from SharedFunctions import *

import cooler
import cooltools

# GLOBAL

JUICER_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../juicer")
GENERATE_SITE_POSITIONS_PATH = os.path.join(JUICER_PATH, "misc/generate_site_positions.py")

# PREPARATION FUNCS

def GetRestrictionSiteLocations(Enzyme, GenomeAssembly, GenomeFA, OutputFile, Logger):
	
	# Logging
	for line in [f"Genome assembly: {GenomeAssembly}", f"Genome FASTA path: {GenomeFA}", f"Enzyme: {Enzyme}", f"Output file: {OutputFile}"]: Logger.info(line)
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = "GetRestrictionSiteLocations",
			Command = f"cd {TempDir}; python3 \"{GENERATE_SITE_POSITIONS_PATH}\" {Enzyme} {GenomeAssembly} \"{GenomeFA}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f" cp {os.path.join(TempDir, '*')} \"{OutputFile}\"",
			Logger = Logger)


# PIPELINE STAGES

def Sra2FastQ(Accession, TopDir, Logger):
	
	# Logging
	for line in [f"Accession: {Accession}", f"Directory: {TopDir}"]: Logger.info(line)
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = "Sra2FastQ",
			Command = f"fastq-dump --split-3 -O {TempDir} \"{Accession}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "GZipFastQ",
			Command = f"gzip {os.path.join(TempDir, '*')}",
			Logger = Logger)
		FastQDir = os.path.join(TopDir, "fastq")
		FastQFiles = glob.glob(os.path.join(TempDir, '*'))
		RenameCommand = '; '.join(["mv " + item + " " + re.sub(r'_([12])\.fastq\.gz$', r'_R\1.fastq.gz', item) for item in FastQFiles if re.match(r".*_([12])\.fastq\.gz$", item) is not None])
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f"mkdir -p \"{FastQDir}\"; {RenameCommand}; cp {os.path.join(TempDir, '*_R*.fastq*')} \"{FastQDir}\"",
			Logger = Logger)


def Juicer(TopDir, Enzyme, RestrictionSiteLocations, GenomeAssembly, GenomeFA, GenomeChromSizes, Threads, Logger):
	
	# Logging
	for line in [f"Directory: {TopDir}", f"Genome assembly: {GenomeAssembly}", f"Genome FASTA path: {GenomeFA}", f"Chrom sizes file path: {GenomeChromSizes}", f"Enzyme: {Enzyme}", f"Restriction sites path: {RestrictionSiteLocations}", f"Threads: {str(Threads)}"]: Logger.info(line)
	
	# Processing
	JuicerScriptPath = os.path.join(JUICER_PATH, "scripts/juicer.sh")
	FullLogText = SimpleSubprocess(
		Name = "Juicer",
		Command = f"bash \"{JuicerScriptPath}\" -t {Threads} -D \"{JUICER_PATH}\" -g {GenomeAssembly} -z \"{GenomeFA}\" -p \"{GenomeChromSizes}\" -s {Enzyme} -y \"{RestrictionSiteLocations}\" -d \"{TopDir}\"",
		Logger = Logger)
	Logger.info(FullLogText.decode('utf-8'))


def MakeInterPairs(InputFile, OutputFile, Logger):
	
	# Logging
	for line in [f"Input file: {InputFile}", f"Output file: {OutputFile}"]: Logger.info(line)
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		TempFile = os.path.join(TempDir, "inter.pairs.txt")
		SimpleSubprocess(
			Name = "MakeInterPairs",
			Command = "awk -F \" \" '{print $2 \"\\t\" $3 \"\\t\" $6 \"\\t\" $7}' \"" + InputFile +"\" > \"" + TempFile + "\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f"cp \"{TempFile}\" \"{OutputFile}\"",
			Logger = Logger)


def Pairs2Cool(PairsFile, CoolFile, GenomeAssembly, GenomeChromSizes, Resolution, Logger):
	
	# Logging
	for line in [f"Input file: {PairsFile}", f"Output file: {CoolFile}", f"Genome assembly: {GenomeAssembly}", f"Chrom sizes file path: {GenomeChromSizes}", f"Resolution [bp]: {Resolution}"]: Logger.info(line)
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		TempFile = os.path.join(TempDir, "inter_no_balanced.cool")
		SimpleSubprocess(
			Name = "MakeCool",
			Command = f"cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 --assembly {GenomeAssembly} \"{GenomeChromSizes}\":{str(Resolution)} \"{PairsFile}\" \"{TempFile}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f"cp \"{TempFile}\" \"{CoolFile}\"",
			Logger = Logger)


def CTaleNormalize(InputFile, OutputFile, Ranges, Logger):
	
	# Logging
	for line in [f"Input file: {InputFile}", f"Output file: {OutputFile}", f"Capture: {Ranges}"]: Logger.info(line)
	
	# Processing
	with tempfile.TemporaryDirectory() as TempDir:
		TempFile = os.path.join(TempDir, "inter.cool")
		SimpleSubprocess(
			Name = "TempFiles",
			Command = f"cp \"{InputFile}\" \"{TempFile}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "NormalizeCool",
			Command = f"ctale_normalize \"{TempFile}\" \"{Ranges}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f"cp \"{TempFile}\" \"{OutputFile}\"",
			Logger = Logger)


def AddWeight(InputFile, Chrom, Resolution, VectorFile, Logger):
	
	# Logging
	for line in [f"Input file: {InputFile}", f"Output file: {VectorFile}", f"Chromosome: {Chrom}", f"Resolution [bp]: {Resolution}"]: Logger.info(line)
	
	# Processing
	try:
		Table = cooler.Cooler(InputFile).bins().fetch(Chrom)
		Weights = Table['weight'].apply(lambda x: 0 if (x == 0) else (1 / x)).fillna(0).rename(f"vector\tc-tale_normalization\t{Chrom}\t{str(Resolution)} BP")
		Weights.to_csv(VectorFile, index=False)
	except Exception as e:
		ErrorMessage = f"Error: {e}"
		Logger.error(ErrorMessage)
		raise RuntimeError(ErrorMessage)
	Logger.info(f"AddWeight successfully finished")


def Vector2HiC(InputFile, OutputFile, VectorFile, Threads, Logger):
	
	# Logging
	for line in [f"Input file: {InputFile}", f"Output file: {OutputFile}", f"Vector: {VectorFile}", f"Threads: {str(Threads)}"]: Logger.info(line)
	
	# Processing
	JuicerToolsPath = os.path.join(JUICER_PATH, "scripts/common/juicer_tools.jar")
	with tempfile.TemporaryDirectory() as TempDir:
		TempFile = os.path.join(TempDir, "inter.hic")
		SimpleSubprocess(
			Name = "TempFiles",
			Command = f"cp \"{InputFile}\" \"{TempFile}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "NormalizeHiC",
			Command = f"java -jar \"{JuicerToolsPath}\" addNorm -j {str(Threads)} \"{TempFile}\" \"{VectorFile}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f"cp \"{TempFile}\" \"{OutputFile}\"",
			Logger = Logger)


# PIPELINE

def BenchmarkPipeline(Name, TopDir, Accessions, Enzyme, RestrictionSiteLocations, GenomeAssembly, GenomeFA, GenomeChromSizes, Capture, Resolution, Threads):
	
	# Make Dirs
	os.mkdir(TopDir)
	TempDir = os.path.join(TopDir, "__temp__")
	os.mkdir(TempDir)
	
	# Logging
	Logger = DefaultLogger(os.path.join(TopDir, "log.txt"))
	Logger.info(f"--- START PROCESSING '{Name}' ---")
	StartTime = time.time()
	
	# Filenames
	Filenames = {
		"MergedNoDups": os.path.join(TempDir, "aligned/merged_nodups.txt"),
		"InterPairs": os.path.join(TempDir, "aligned/inter.pairs.txt"),
		"InterNoBalanced": os.path.join(TempDir, "aligned/inter_no_balanced.cool"),
		"InterCool": os.path.join(TempDir, "aligned/inter.cool"),
		"Vector": os.path.join(TempDir, "aligned/vector.txt"),
		"InterHic": os.path.join(TempDir, "aligned/inter.hic"),
		"InterHicNormalized": os.path.join(TempDir, "aligned/inter_norm.hic"),
		"InterStat": os.path.join(TempDir, "aligned/inter.txt")
		}
	
	CopyFilenames = {
		Filenames["InterCool"]: os.path.join(TopDir, "inter.cool"),
		Filenames["InterHicNormalized"]: os.path.join(TopDir, "inter.hic"),
		Filenames["InterPairs"]: os.path.join(TopDir, "inter.pairs.txt"),
		Filenames["InterNoBalanced"]: os.path.join(TopDir, "inter_no_balanced.cool"),
		Filenames["InterStat"]: os.path.join(TopDir, "inter_statistics.txt"),
		Filenames["MergedNoDups"]: os.path.join(TopDir, "merged_nodups.txt"),
		Filenames["Vector"]: os.path.join(TopDir, "vector.txt")
		}
	
	# Pipeline
	for Accession in Accessions:
		Sra2FastQ(
			Accession = Accession,
			TopDir = TempDir,
			Logger = Logger)
	
	Juicer(
		TopDir = TempDir,
		Enzyme = Enzyme,
		RestrictionSiteLocations = RestrictionSiteLocations,
		GenomeAssembly = GenomeAssembly,
		GenomeFA = GenomeFA,
		GenomeChromSizes = GenomeChromSizes,
		Threads = Threads,
		Logger = Logger)
	
	MakeInterPairs(
		InputFile = Filenames["MergedNoDups"],
		OutputFile = Filenames["InterPairs"],
		Logger = Logger)
	
	Pairs2Cool(
		PairsFile = Filenames["InterPairs"],
		CoolFile = Filenames["InterNoBalanced"],
		GenomeAssembly = GenomeAssembly,
		GenomeChromSizes = GenomeChromSizes,
		Resolution = Resolution,
		Logger = Logger)
	
	CTaleNormalize(
		InputFile = Filenames["InterNoBalanced"],
		OutputFile = Filenames["InterCool"],
		Ranges = f"{Capture[0]}:{Capture[1]:,}-{Capture[2]:,}",
		Logger = Logger)
	
	AddWeight(
		InputFile = Filenames["InterCool"],
		Chrom = Capture[0],
		Resolution = Resolution,
		VectorFile = Filenames["Vector"],
		Logger = Logger)
	
	Vector2HiC(
		InputFile = Filenames["InterHic"],
		OutputFile = Filenames["InterHicNormalized"],
		VectorFile = Filenames["Vector"],
		Threads = Threads,
		Logger = Logger)
	
	for source, dest in CopyFilenames.items():
		SimpleSubprocess(
			Name = "CopyResult",
			Command = f"cp \"{source}\" \"{dest}\"",
			Logger = Logger)
	
	Logger.info(f"'{Name}' FINISHED, SUMMARY TIME - %s" % (SecToTime(time.time() - StartTime)))

# --------------------------

Data = pandas.read_csv("/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Project/Blocks.tsv", sep='\t', comment='#')
for index, line in Data.iterrows():
	Name = line["Accession"].replace(";", "-")
	TopDir = f"/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Project/{Name}"
	Accessions = [f"/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Source/{item}.sra" for item in line["Accession"].split(';')]
	GenomeFA = f"/Data/DataBases/{line['Genome']}/{line['Genome']}_canonic.fa"
	GenomeChromSizes = f"/Data/DataBases/{line['Genome']}/{line['Genome']}_canonic.chrom.sizes"
	RestrictionSiteLocations = f"/Data/UserData/FairWind/Ya.Cloud/core/pipeline/data/restriction_sites/RestrictionSites_{line['Enzyme']}_{line['Genome']}.txt"
	Capture = [line["Chrom"], line["Start"], line["End"]]
	Resolution = 5000
	Threads = 12
	BenchmarkPipeline(
		Name = Name,
		TopDir = TopDir,
		Accessions = Accessions,
		Enzyme = line['Enzyme'],
		RestrictionSiteLocations = RestrictionSiteLocations,
		GenomeAssembly = line['Genome'],
		GenomeFA = GenomeFA,
		GenomeChromSizes = GenomeChromSizes,
		Capture = Capture,
		Resolution = Resolution,
		Threads = Threads)
