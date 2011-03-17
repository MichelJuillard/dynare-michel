if HAVE_PDFLATEX
if HAVE_BIBTEX
pdf-local: gsa.pdf
endif
endif

SRC = gsa.tex marco.bib

EXTRA_DIST = $(SRC)

gsa.pdf: $(SRC)
	$(PDFLATEX) gsa
	$(BIBTEX) gsa
	$(PDFLATEX) gsa
	$(PDFLATEX) gsa

clean-local:
	rm -f *.pdf *.log *.aux *.toc *.lof *.blg *.bbl *.out *~
