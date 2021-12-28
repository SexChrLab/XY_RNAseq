snakemake --snakefile gencode.smk --rulegraph -n -p all > rulegraph.txt
cat rulegraph.txt | dot -Tsvg > rulegraph.svg
git add gencode.smk rulegraph.txt rulegraph.svg make_rulegraph.sh


snakemake --snakefile gencode.smk --dag -n -p all > dag.txt
cat dag.txt | dot -Tsvg > dag.svg
git add gencode.smk dag.txt dag.svg make_dag.sh

snakemake --snakefile gencode.smk --filegraph -n -p all > filegraph.txt
cat filegraph.txt | dot -Tsvg > filegraph.svg

git add gencode.smk filegraph.txt filegraph.svg make_rulegraph.sh
git commit -m "rulegraph filegraph dag creation"
