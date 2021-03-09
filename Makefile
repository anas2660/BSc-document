
pdf: BSc.tex
	pdflatex BSc.tex

doc: pdf
	cp -f BSc.pdf Doc.pdf

clean:
	rm -f BSc.aux BSc.toc BSc.log
