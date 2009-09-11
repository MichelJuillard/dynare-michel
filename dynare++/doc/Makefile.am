EXTRA_DIST = dynare++-ramsey.tex dynare++-tutorial.tex

if HAVE_PDFLATEX
pdf-local: dynare++-ramsey.pdf dynare++-tutorial.pdf
endif

%.pdf: %.tex
	$(PDFLATEX) $<
	$(PDFLATEX) $<
	$(PDFLATEX) $<

CLEANFILES = *.pdf *.log *.aux *.out *.toc


