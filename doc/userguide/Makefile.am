if HAVE_PDFLATEX
if HAVE_BIBTEX
all-local: UserGuide.pdf
endif
endif

SRC = UserGuide.tex Graphics/DynareTitle.pdf DynareBib.bib \
	ch-intro.tex ch-inst.tex ch-solbase.tex ch-soladv.tex ch-estbase.tex \
	ch-estadv.tex ch-solbeh.tex ch-estbeh.tex ch-ramsey.tex ch-trouble.tex

EXTRA_DIST = $(SRC)

UserGuide.pdf: $(SRC)
	$(PDFLATEX) UserGuide
	$(BIBTEX) UserGuide
	$(PDFLATEX) UserGuide
	$(PDFLATEX) UserGuide

clean-local:
	rm -f *.log *.aux *.toc *.lof *.blg *.bbl *.out
	rm -f UserGuide.pdf
