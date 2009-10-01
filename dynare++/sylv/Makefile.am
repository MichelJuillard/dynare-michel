SUBDIRS = cc testing

EXTRA_DIST = sylvester.tex change_log.html matlab

if HAVE_PDFTEX
pdf-local: sylvester.pdf

sylvester.pdf: sylvester.tex
	$(PDFTEX) sylvester
endif

CLEANFILES = *.pdf *.log
