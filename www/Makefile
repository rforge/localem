REXE = R
PANDOC = pandoc

all: simulations.html

%.html: %.md
	$(PANDOC) --to=html5 --output=$@ $<

%.md: %.Rmd
	$(REXE) -e "knitr::knit('$<', encoding='UTF-8')" highRes
	
simulations.Rmd: FORCE
	$(REXE) -e "file.copy(system.file(file.path('doc','kentucky.Rmd'),package='localEM'),'simulations.Rmd')"

FORCE:	

clean:
	rm *.md *.Rmd

