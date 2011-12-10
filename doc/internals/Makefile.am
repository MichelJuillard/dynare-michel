EXTRA_DIST = dynare-internals.org

if ENABLE_ORG_EXPORT
html-local:
	emacs --batch --visit=dynare-internals.org --funcall org-export-as-html-batch
endif

clean-local:
	rm -rf *.html ltxpng
