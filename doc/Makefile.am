SUBDIRS = preprocessor macroprocessor userguide parallel internals gsa

info_TEXINFOS = dynare.texi

if HAVE_TEXI2HTML
if HAVE_LATEX2HTML
html-local: dynare.texi
	rm -rf dynare.html
	mkdir -p dynare.html
	cd dynare.html && $(TEXI2HTML) --l2h --split section --prefix index ../dynare.texi
endif
endif

PDF_TARGETS =

if HAVE_PDFLATEX
PDF_TARGETS += guide.pdf bvar-a-la-sims.pdf
if HAVE_BIBTEX
PDF_TARGETS += dr.pdf
endif
endif

pdf-local: $(PDF_TARGETS)

EXTRA_DIST = guide.tex guide.bbl bibmad.sty bvar-a-la-sims.tex dr.tex dr.bib

guide.pdf: guide.tex guide.bbl bibmad.sty
	$(PDFLATEX) guide
	$(PDFLATEX) guide

bvar-a-la-sims.pdf: bvar-a-la-sims.tex
	$(PDFLATEX) bvar-a-la-sims
	$(PDFLATEX) bvar-a-la-sims

dr.pdf: dr.tex
	$(PDFLATEX) dr
	$(BIBTEX) dr
	$(PDFLATEX) dr
	$(PDFLATEX) dr

clean-local:
	rm -f *~ *.pdf *.log *.aux *.out *.blg *.bbl
