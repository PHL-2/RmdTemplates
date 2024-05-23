# output gzipped file
zcat "$1" | 
#get every 2nd line in fastq file and include the sample name (modified filename)
  awk -v filename="$1" 'NR % 4 == 2 {gsub(".*scrubbed\\.[0-9]+\\.|_S.*", "" ,filename); print filename" "$1}' |
#get top nth most common reads
  sort | uniq -c | sort -r | head -n $2 |\
#make into fasta file for BLAST upload
  cat -n | sed 's/^[\t ]*/>/g' | sed 's/[\t ]\+/_/g' | sed 's/\(.*\)_/\1\n/g' > "$3"