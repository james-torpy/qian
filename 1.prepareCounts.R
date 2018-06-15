### 10.prepareCounts.R ###

# breaks gcCounts up into dataframe for each class and
# saves each as RData objects #

# define starting variables:
project <- "qian"
expName <- "exp1"

# define directories:
homeDir <- "/share/ScratchGeneral/jamtor/"
#homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", 
  project, "/")
resultsDir <- paste0(projectDir, "/results/", expName, "/")
RobjectDir <- paste0(projectDir, "/Robjects/", 
  expName, "/")
inDir <- paste0(resultsDir, "/htseq/")
rawDir <- paste0(projectDir, 
  "/raw_files/", expName, "/fastq/")
starGC_dir <- paste0(resultsDir, "/star/GC/")
starRibo_dir <- paste0(resultsDir, "/star/ribo/", 
  expName)

system(paste0("mkdir -p ", RobjectDir))



### 1. Check for and discard bad inputs ###

# write function to ensure all htseq output files have same number of lines:
file_check <- function(type) {
  in_files <- list.files(inDir, pattern = type, full.names = T)
  good_files <- c()
  bad_files <- c()
  for ( i in 1:length(in_files) ) {
    if (i==1) {
      line_no <- as.numeric(
        strsplit(
          system(
            paste0("wc -l ", in_files[i]), intern=T
          ), " "
        )[[1]][1]
      )
      print(paste0(in_files[i], " has ", line_no, " lines"))
      good_files <- append(in_files[i], good_files)
    } else {
      current_no <- as.numeric(
        strsplit(
          system(
            paste0("wc -l ", in_files[i]), intern=T
          ), " "
        )[[1]][1]
      )
      if (current_no == line_no) {
        print(paste0(in_files[i], " has the same number of lines as first htseq output file ",
                     in_files[1]))
        good_files <- append(in_files[i], good_files)
      } else {
        print(paste0(in_files[i], " has a different number of lines as first htseq output file ",
                     in_files[1]))
        bad_files <- append(in_files[i], bad_files)
      }
    }
  }
  res <- list(good_files, bad_files)
  names(res) <- c("good_files", "bad_files")
  
  return(res)
}

gc_file_check <- file_check("gc")
gc_files <- gc_file_check$good_files
gc_IDs <- gsub(
  "\\..*$", "", basename(gc_files)
)

all_file_check <- file_check("all")
all_files <- all_file_check$good_files
all_IDs <- gsub(
  "\\..*$", "", basename(all_files)
)

all_reps_file_check <- file_check("all")
all_reps_files <- all_reps_file_check$good_files
all_reps_IDs <- gsub(
  "\\..*$", "", basename(all_reps_files)
)

# fetch those counts file ids which worked for both all and gc:
IDs <- Reduce(intersect, list(gc_IDs, all_IDs, all_reps_IDs))


### 2. Load in inputs and format ###

# create function to load inputs and format:
load_and_format <- function(id, type) {
  for ( j in 1:length(id) ) {
    print(paste0("loading ", id[j]))
    if (j==1) {
      counts <- data.frame(read.table(file=paste0(inDir, "/", IDs[j], ".", type, ".htseq.txt")))
    } else {
      counts <- cbind(counts, data.frame(read.table(file=paste0(inDir, "/", IDs[j], ".", type, ".htseq.txt")))[,2])
    }
  }
  colnames(counts) <- c("gene_id", IDs)
  saveRDS(counts, paste0(RobjectDir, type, "_counts_temp.rds"))
  
  if ( type == "gc") {
    # aggregate multiple types of same gene:
    counts$gene_id <- gsub("\\..*$", "", counts$gene_id)
    counts <- aggregate(.~gene_id, counts, mean)
  }
  # remove duplicate rows:
  counts <- unique(counts)
  # remove specs lines:
  counts <- counts[grep("__", counts$gene_id, invert=T),]

  if ( type == "all" ) {
    counts <- apply(counts[,2:ncol(counts)], 2, sum)
  }
  
  return(counts)
}

# load and format gc counts:
gc_counts <- load_and_format(IDs, "gc")
# save the counts:
if (!file.exists(paste0(RobjectDir, "/gc_counts.htseq.rds"))) {
  saveRDS(gc_counts, file=paste0(RobjectDir, "/gc_counts.htseq.rds"))
}

# load and format all counts:
all_counts <- load_and_format(IDs, "all")
# save the counts:
if (!file.exists(paste0(RobjectDir, "/all_counts.htseq.rds"))) {
  saveRDS(all_counts, file=paste0(RobjectDir, "/all_counts.htseq.rds"))
}

# load and format all repeat counts:
all_rep_counts <- load_and_format(IDs, "all_reps")
# save the counts:
if (!file.exists(paste0(RobjectDir, "/all_rep_counts.htseq.rds"))) {
  saveRDS(all_rep_counts, file=paste0(RobjectDir, "/all_rep_counts.htseq.rds"))
}
