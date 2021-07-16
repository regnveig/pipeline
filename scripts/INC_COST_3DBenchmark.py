from SharedFunctions import *

JUICER_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../juicer")
GENERATE_SITE_POSITIONS_PATH = os.path.join(JUICER_PATH, "misc/generate_site_positions.py")

def Sra2FastQ(Accession, TopDir, Logger):
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

def GetRestrictionSiteLocations(Enzyme, GenomeAssembly, GenomeFA, OutputFile, Logger):
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = "GetRestrictionSiteLocations",
			Command = f"cd {TempDir}; python3 \"{GENERATE_SITE_POSITIONS_PATH}\" {Enzyme} {GenomeAssembly} \"{GenomeFA}\"",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f" cp {os.path.join(TempDir, '*')} \"{OutputFile}\"",
			Logger = Logger)

def Juicer(TopDir, Enzyme, RestrictionSiteLocations, GenomeAssembly, GenomeFA, GenomeChromSizes, Threads):
	JuicerScriptPath = os.path.join(JUICER_PATH, 'scripts/juicer.sh')
	FullLogText = SimpleSubprocess(
		Name = "Juicer",
		Command = f"bash \"{JuicerScriptPath}\" -t {Threads} -D \"{JUICER_PATH}\" -g {GenomeAssembly} -z \"{GenomeFA}\" -p \"{GenomeChromSizes}\" -s {Enzyme} -y \"{RestrictionSiteLocations}\" -d \"{TopDir}\"",
		Logger = Logger)
	Logger.info(FullLogText)

def MakeInterPairs(InputFile, OutputFile):
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

Logger = DefaultLogger("/dev/null")
#Sra2FastQ("/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Source/SRR5453505.sra", "/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Project/TEST", Logger)
#Juicer("/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Project/TEST", "DpnII", "/Data/UserData/FairWind/Ya.Cloud/core/pipeline/data/restriction_sites/RestrictionSites_DpnII_mm10.txt", "mm10", "/Data/DataBases/mm10/mm10_canonic.fa", "/Data/DataBases/mm10/mm10_canonic.chrom.sizes", 12)
MakeInterPairs("/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Project/TEST/aligned/merged_nodups.txt", "/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Project/TEST/aligned/inter.pairs.txt")
#def t(x): GetRestrictionSiteLocations(x[0], x[1], f"/Data/DataBases/{x[1]}/{x[1]}_canonic.fa", f"/Data/UserData/FairWind/Ya.Cloud/core/pipeline/data/restriction_sites/RestrictionSites_{x[0]}_{x[1]}.txt", Logger)
#with Threading("Restrict", Logger, 4) as pool: pool.map(t, [["DpnII", "mm10"], ["DpnII", "hg38"], ["HindIII", "mm10"], ["HindIII", "hg38"]])
#for Site in ["DpnII", "HindIII"]:
	#for genome in ["mm10", "hg38"]:
		#GetRestrictionSiteLocations(Site, genome, f"/Data/DataBases/{genome}/{genome}_canonic.fa", f"/Data/UserData/FairWind/Ya.Cloud/core/pipeline/data/restriction_sites/RestrictionSites_{Site}_{genome}.txt", Logger)
