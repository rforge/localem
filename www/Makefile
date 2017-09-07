REXE = R
PANDOC = pandoc

all: kentucky.html kentuckyHighRes.html

.PRECIOUS: *.md

kentucky.Rmd: ../pkg/localem/inst/doc/kentucky.Rmd
	cp $< $@
	
kentuckyHighRes.md: kentucky.Rmd
	$(REXE) -e "knitr::knit('$<', '$@', encoding='UTF-8')" highRes

kentucky.md: kentucky.Rmd
	$(REXE) -e "knitr::knit('$<', '$@', encoding='UTF-8')"

%.html: %.md
	$(PANDOC) --standalone --smart --filter=pandoc-citeproc --filter=pandoc-crossref --from=markdown --to=html5 --output=$@ $<

clean:
	rm *.aux *.blg *.log *.nav *.out *.snm *.toc *.tex *.vrb *.gri *.grd
	
