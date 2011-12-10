if HAVE_PDFLATEX
if HAVE_BIBTEX
pdf-local: parallel.pdf
endif
endif

SRC = AvenueParadigm.pdf iVaNo_gain.pdf iVaNo_time_comp.pdf marco.bib					\
	netbook_complete_comp.pdf netbook_complete_openclose.pdf										\
	netbook_partial_comp.pdf netbook_partial_openclose.pdf parallel.tex					\
	quest_complete_comp.pdf quest_complete_openclose.pdf quest_partial_comp.pdf	\
	quest_partial_openclose.pdf RWMH_quest1_PriorsAndPosteriors1Comp.pdf				\
	RWMH_quest1_PriorsAndPosteriors2Comp.pdf																		\
	RWMH_quest1_PriorsAndPosteriors3Comp.pdf																		\
	RWMH_quest1_PriorsAndPosteriors4Comp.pdf																		\
	RWMH_quest1_PriorsAndPosteriors5Comp.pdf																		\
	RWMH_quest1_PriorsAndPosteriors6Comp.pdf																		\
	RWMH_quest1_PriorsAndPosteriors7Comp.pdf waitbars1.pdf waitbars2.pdf				\
	waitbarsP.pdf

EXTRA_DIST = $(SRC)

parallel.pdf: $(SRC)
	$(PDFLATEX) parallel
	$(BIBTEX) parallel
	$(PDFLATEX) parallel
	$(PDFLATEX) parallel
	$(PDFLATEX) parallel

clean-local:
	rm -f *.log *.aux *.toc *.blg *.bbl *.out
	rm -f parallel.pdf
