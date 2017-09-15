
REXE = R
PANDOC = pandoc

all: kentucky.html

.PRECIOUS: *.md

%.html: %.md
	$(PANDOC) --standalone --smart --filter=pandoc-citeproc --filter=pandoc-crossref --from=markdown --to=html5 --output=$@ $<

kentuckyHighRes.md: kentucky.Rmd
	$(REXE) -e "knitr::knit('$<', '$@', encoding='UTF-8')" highRes

kentucky.md: kentucky.Rmd
	$(REXE) -e "knitr::knit('$<', '$@', encoding='UTF-8')"

clean:
	rm *.aux *.blg *.log *.nav *.out *.snm *.toc *.tex *.vrb *.gri *.grd
	
