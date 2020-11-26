# Output accession numbers from a query file which are not present in a fasta file
# USAGE: awk -f find_missing_accessions.awk accessions_file fasta_file
# 
# Accessions_file should have one accession per line
# If all accessions are found in the fasta file, there will be no output

BEGIN { 
    FS=" "                               # break lines into separate words
    unique_count=1
}

(NR==FNR) {                              # file 1 (accessions_file)
    if (!($1 in seen)) {
        seen[$1]++                       # lines, indeced by string
        lines1[unique_count]=$1          # lines, indexed by order
        unique_count++
    }
}

(NR!=FNR) {                              # file 2 (fasta_file)
    acc=substr($1, 2)                    # taking first word on each line, minus the first char (e.g. "abc" in ">abc dce ...")
    lines2[acc]++                        # indexed by string
}

END {
    header=0
    
    for (i=1; i<=unique_count; i++) {
        if (!(lines1[i] in lines2) && lines1[i] != "") {
            if (header==0) {
                print "Accessions not in fasta file:"
                header=1
            } 
            
            print lines1[i]
        }
    }
}