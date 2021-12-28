snakemake --snakefile ensembl.smk --rulegraph -n -p all > rulegraph.txt
cat rulegraph.txt | dot -Tsvg > rulegraph.svg
git add ensembl.smk rulegraph.txt rulegraph.svg make_rulegraph.sh
git commit -m "rulegraph creation"
