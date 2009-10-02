SUBDIRS = cc testing

EXTRA_DIST = sylvester.tex change_log.html matlab

if HAVE_PDFETEX
pdf-local: sylvester.pdf

sylvester.pdf: sylvester.tex
	$(PDFETEX) sylvester
endif

CLEANFILES = *.pdf *.log
