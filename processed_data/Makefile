all: feature-table.tsv

%.tsv: %.biom
	biom convert -i $< --to-tsv -o $@
