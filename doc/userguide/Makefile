all: UserGuide.pdf

SRC = UserGuide.tex Graphics/DynareTitle.pdf DynareBib.bib \
	ch-intro.tex ch-inst.tex ch-solbase.tex ch-soladv.tex ch-estbase.tex \
	ch-estadv.tex ch-solbeh.tex ch-estbeh.tex ch-ramsey.tex ch-trouble.tex

UserGuide.pdf: 
	pdflatex UserGuide
	bibtex UserGuide
	pdflatex UserGuide
	pdflatex UserGuide

clean:
	rm -f *.log *.aux *.toc *.lof *.blg *.bbl *.out
	rm -f UserGuide.pdf
