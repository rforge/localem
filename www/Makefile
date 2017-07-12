REXE = R
PANDOC = pandoc

all: simulations.html

.PRECIOUS: *.md


%.html: %.md
	$(PANDOC) --standalone --smart --filter=pandoc-citeproc --filter=pandoc-crossref --from=markdown --to=html5 --output=$@ $<

%.md: %.Rmd
	$(REXE) -e "knitr::knit('$<', encoding='UTF-8')" highRes
	
simulations.Rmd:
	$(REXE) -e "file.copy(system.file(file.path('doc','kentucky.Rmd'),package='localEM'),'simulations.Rmd')"

clean:
	rm *.md *.Rmd

