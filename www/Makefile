REXE = R
PANDOC = pandoc
LATEX = xelatex
SED = sed

all: kentucky.html kentuckyHighRes.html

toronto: toronto.html

.PRECIOUS: *.md

%.Rmd: ../pkg/localem/inst/extdoc/%.Rmd
	cp $< $@
	
kentuckyHighRes.md: kentucky.Rmd
	$(REXE) -e "knitr::knit('$<', '$@', encoding='UTF-8')" highRes

%.md: %.Rmd
	$(REXE) -e "knitr::knit('$<', '$@', encoding='UTF-8')"

%.html: %.md
	$(PANDOC) --standalone --filter=pandoc-citeproc --filter=pandoc-crossref --from=markdown --to=html5 --output=$@ $<

%.tex: %.md
	$(PANDOC) --standalone --biblatex --filter=pandoc-citeproc --filter=pandoc-crossref --from=markdown --to=latex $< | $(SED) s/\\\\usepackage{subfig}// > $@

%.pdf: %.tex
	$(LATEX) $<

clean:
	rm *.aux *.blg *.log *.nav *.out *.snm *.toc *.tex *.vrb *.gri *.grd
	
