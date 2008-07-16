all: guide.pdf manual-html/index.html manual.pdf bvar-a-la-sims.pdf
	make -C preprocessor
	make -C macroprocessor

guide.pdf: guide.tex guide.bbl bibmad.sty
	pdflatex guide
	pdflatex guide

bvar-a-la-sims.pdf: bvar-a-la-sims.tex
	pdflatex bvar-a-la-sims
	pdflatex bvar-a-la-sims

manual-html/index.html: manual.xml dynare_html.xsl
	# Make sure that you have installed the docbook-xsl package, otherwise xsltproc will take a very long time
	xsltproc -stringparam base.dir ./manual-html/ dynare_html.xsl manual.xml 

manual.pdf: manual.xml
	dblatex manual.xml

clean:
	rm -f *~
	rm -f guide.{pdf,log,aux}
	rm -f bvar-a-la-sims.{pdf,log,aux}
	rm -rf manual.pdf manual-html
	make -C preprocessor clean
	make -C macroprocessor clean
