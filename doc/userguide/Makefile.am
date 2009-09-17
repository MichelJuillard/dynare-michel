if HAVE_PDFLATEX
if HAVE_BIBTEX
pdf-local: UserGuide.pdf
endif
endif

SRC = UserGuide.tex Graphics/DynareTitle.pdf DynareBib.bib \
	ch-intro.tex ch-inst.tex ch-solbase.tex ch-soladv.tex ch-estbase.tex \
	ch-estadv.tex ch-solbeh.tex ch-estbeh.tex ch-ramsey.tex ch-trouble.tex \
	P_DynareStruct2.pdf P_flowest.pdf P_MH2.pdf P_ModStruct2.pdf P_ModStruct3.pdf \
	P_ModStruct4.pdf P_ModStruct5.pdf P_SchorfMod.pdf P_ShockModel2.pdf

EXTRA_DIST = $(SRC) Graphics models

UserGuide.pdf: $(SRC)
	$(PDFLATEX) UserGuide
	$(BIBTEX) UserGuide
	$(PDFLATEX) UserGuide
	$(PDFLATEX) UserGuide
	$(PDFLATEX) UserGuide

clean-local:
	rm -f *.log *.aux *.toc *.lof *.blg *.bbl *.out
	rm -f UserGuide.pdf
