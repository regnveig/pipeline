from SharedFunctions import *

def Sra2FastQ(Accession, OutputFolder, Logger):
	with tempfile.TemporaryDirectory() as TempDir:
		SimpleSubprocess(
			Name = "Sra2FastQ",
			Command = f"fastq-dump --split-3 -O {TempDir} {Accession}",
			Logger = Logger)
		SimpleSubprocess(
			Name = "GZipFastQ",
			Command = f"gzip {os.path.join(TempDir, '*')}",
			Logger = Logger)
		SimpleSubprocess(
			Name = "CopyFiles",
			Command = f"cp -R {os.path.join(TempDir, '*')} \"{OutputFolder}\"",
			Logger = Logger)


Logger = DefaultLogger("/dev/null")
Sra2FastQ("/Data/NGS_Data/20210714_INC_COST_3DBenchmark/Source/SRR5453505.sra", "./", Logger)
