
pdf: BSc.pdf

BSc.pdf: BSc.tex BSc.bbl
	lualatex BSc.tex -interaction=nonstopmode -halt-on-error

doc: pdf
	cp -f BSc.pdf Doc.pdf

BSc.bbl: References.bib
	biber BSc

clean:
	rm -rf BSc.aux BSc.toc BSc.log BSc.run.xml BSc.bcf _region_* BSc.blg BSc.bbl prv_BSc.log
