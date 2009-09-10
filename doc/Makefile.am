SUBDIRS = preprocessor macroprocessor userguide

TARGETS =

if HAVE_PDFLATEX
TARGETS += guide.pdf bvar-a-la-sims.pdf
endif

if HAVE_BIBTEX
TARGETS += manual.pdf
endif

if HAVE_XSLT
TARGETS += manual-html/index.html
endif

all-local: $(TARGETS)

EXTRA_DIST = guide.tex guide.bbl bibmad.sty bvar-a-la-sims.tex manual.xml dynare_html.xsl manual.xml

guide.pdf: guide.tex guide.bbl bibmad.sty
	$(PDFLATEX) guide
	$(PDFLATEX) guide

bvar-a-la-sims.pdf: bvar-a-la-sims.tex
	$(PDFLATEX) bvar-a-la-sims
	$(PDFLATEX) bvar-a-la-sims

manual-html/index.html: manual.xml dynare_html.xsl
	$(XSLTPROC) -stringparam base.dir ./manual-html/ dynare_html.xsl manual.xml

manual.pdf: manual.xml
	$(DBLATEX) manual.xml

clean-local:
	rm -f *~ *.pdf *.log *.aux
	rm -rf manual-html
