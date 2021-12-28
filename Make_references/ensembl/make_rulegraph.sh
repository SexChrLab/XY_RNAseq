snakemake --snakefile ensembl.smk --rulegraph -n -p all > rulegraph.txt
cat rulegraph.txt | dot -Tsvg > rulegraph.svg
git add ensembl.smk rulegraph.txt rulegraph.svg make_rulegraph.sh
git commit -m "rulegraph creation"

snakemake --snakefile ensembl.smk --dag -n -p all > dag.txt
cat dag.txt | dot -Tsvg > dag.svg
git add ensembl.smk dag.txt dag.svg make_rulegraph.sh

snakemake --snakefile ensembl.smk --filegraph -n -p all > filegraph.txt
cat filegraph.txt | dot -Tsvg > filegraph.svg
git add ensembl.smk filegraph.txt filegraph.svg make_rulegraph.sh

git commit -m "rulegraph filegraph dag creation"
