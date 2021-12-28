snakemake --snakefile gencode.smk --rulegraph -n -p all > rulegraph.txt
cat rulegraph.txt | dot -Tsvg > rulegraph.svg
git add gencode.smk rulegraph.txt rulegraph.svg make_rulegraph.sh
git commit -m "rulegraph creation"
