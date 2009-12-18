SUBDIRS = preprocessor macroprocessor userguide

PDF_TARGETS =

if HAVE_PDFLATEX
PDF_TARGETS += guide.pdf bvar-a-la-sims.pdf dr.pdf
endif

if HAVE_DBLATEX
PDF_TARGETS += manual.pdf
endif

pdf-local: $(PDF_TARGETS)

if HAVE_XSLT
html-local: manual-html/index.html
endif

EXTRA_DIST = guide.tex guide.bbl bibmad.sty bvar-a-la-sims.tex dr.tex manual.xml dynare_html.xsl dynare.xsl

guide.pdf: guide.tex guide.bbl bibmad.sty
	$(PDFLATEX) guide
	$(PDFLATEX) guide

bvar-a-la-sims.pdf: bvar-a-la-sims.tex
	$(PDFLATEX) bvar-a-la-sims
	$(PDFLATEX) bvar-a-la-sims

dr.pdf: dr.tex
	$(PDFLATEX) dr
	$(PDFLATEX) dr

manual-html/index.html: manual.xml dynare_html.xsl dynare.xsl
	$(XSLTPROC) -stringparam base.dir ./manual-html/ dynare_html.xsl manual.xml

manual.pdf: manual.xml
	$(DBLATEX) manual.xml

clean-local:
	rm -f *~ *.pdf *.log *.aux *.out
	rm -rf manual-html
