all: guide.pdf manual/index.html #manual.pdf

guide.pdf: guide.tex
	pdflatex guide

manual/index.html: manual.xml
	xsltproc -stringparam base.dir ./manual/ dynare_html.xsl manual.xml 


manual.pdf: manual.xml
	xsltproc dynare_fo.xsl manual.xml >manual.fo
	pdfxmltex manual.fo	
	pdfxmltex manual.fo	

install: manual.pdf manual/index.html
	scp manual.pdf pythie.cepremap.cnrs.fr:public_html/mambo/download/manual/
	scp manual/* pythie.cepremap.cnrs.fr:public_html/mambo/download/manual/
