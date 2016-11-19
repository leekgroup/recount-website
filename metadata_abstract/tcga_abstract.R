## Text taken from https://cancergenome.nih.gov/abouttcga

abstracts <- data.frame(
    study_abstract = "The Cancer Genome Atlas (TCGA) is a collaboration between the National Cancer Institute (NCI) and the National Human Genome Research Institute (NHGRI) that has generated comprehensive, multi-dimensional maps of the key genomic changes in 33 types of cancer. The TCGA dataset, comprising more than two petabytes of genomic data, has been made publically available, and this genomic information helps the cancer research community to improve the prevention, diagnosis, and treatment of cancer.",
    study_accession = 'TCGA',
    stringsAsFactors = FALSE)

save(abstracts, file = 'abstracts_tcga.Rdata')
